###########################################################Argparse###########################################################

library(argparse)

# Add and define the --model, aa_genes, --survival_data, and --outdir arguments
parser <- ArgumentParser()

parser$add_argument('--model', required=TRUE, help="Path to the model matrix CSV. 
                    Required cols: sample, ecdna, subtype, group, sample_type, method.")
parser$add_argument('--aa_genes', required=TRUE, help="Path to the AmpliconArchitect gene summary CSV.")
parser$add_argument('--survival_data', required=TRUE, help="Path to the survival_data CSV.
                    Required cols: sample, time, censored.")
parser$add_argument('--outdir', required=TRUE, help="Path to output directory.")
args = parser$parse_args()
model = args$model
aa_genes = args$aa_genes
survival_data = args$survival_data
outdir = args$outdir

###############################################################################################################################

###########################################################Load Libraries###########################################################

library(GenomicFeatures)
library(org.Hs.eg.db)
library(immunedeconv)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(maftools)
library(msigdbr)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(minfi)
library(survival)
library(survminer)
library(DESeq2)
library(gcap)
library(clusterProfiler)

# Source functions for use in main.R
source("/Users/drainford/Desktop/skcm_final/container/functions.R")

###############################################################################################################################

###################################################### Group Definitions ######################################################
# Read in model matrix from the provided file path
model_matrix = read.csv(model)

# Define experimental groups based on subtype and ecDNA status
# Subtypes: TWT (Triple Wild Type), BNN (BRAF, NRAS, NF1)
twt_samples = model_matrix[model_matrix$subtype == "TWT", ]$sample  # Samples for TWT subtype
bnn_samples = model_matrix[model_matrix$subtype == "BNN", ]$sample  # Samples for BNN subtype
ecdna_pos = model_matrix[model_matrix$ecdna == "Positive", ]$sample  # ecDNA positive samples
ecdna_neg = model_matrix[model_matrix$ecdna == "Negative", ]$sample  # ecDNA negative samples

# Define groups based on both subtype and ecDNA status
bnn_neg = model_matrix[model_matrix$group == "BNN_ecDNA-Negative", ]$sample  # BNN ecDNA-Negative samples
bnn_pos = model_matrix[model_matrix$group == "BNN_ecDNA-Positive", ]$sample  # BNN ecDNA-Positive samples
twt_neg = model_matrix[model_matrix$group == "TWT_ecDNA-Negative", ]$sample  # TWT ecDNA-Negative samples
twt_pos = model_matrix[model_matrix$group == "TWT_ecDNA-Positive", ]$sample  # TWT ecDNA-Positive samples

# Define groups by method
aa_samples = model_matrix[model_matrix$method == "AA", ]$sample  # Amplicon Architect (AA) samples
gcap_samples = model_matrix[model_matrix$method == "GCAP", ]$sample  # GCAP samples
#############################################################################################################################


####################################################### ecDNA Frequency #######################################################

# Calculate the frequency of ecDNA by subtype
ecdna_counts = model_matrix %>%
  group_by(subtype, ecdna) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate ecDNA frequency for BNN and TWT subtypes
# Frequency is calculated as the ratio of positive ecDNA samples over the total number of samples in each subtype
bnn = as.numeric(round(ecdna_counts[2,3] / (ecdna_counts[2,3] + ecdna_counts[1,3]),2) * 100)
twt = as.numeric(round(ecdna_counts[4,3] / (ecdna_counts[4,3] + ecdna_counts[3,3]),2) * 100)
ecdna_freq = data.frame(sample = c("BNN", "TWT"), freq = c(bnn, twt))

write.csv(ecdna_freq, paste0(outdir, "ecdna_frequency.csv"))

#############################################################################################################################

#################################################### Metastatic Frequency #####################################################

# Filter out unknown sample types to retain only Primary Tumor and Metastatic types
meta_matrix = model_matrix[model_matrix$sample_type %in% c("Primary Tumor", "Metastatic"), ]

# Get counts of metastatic and primary tumors for each group
sample_counts = meta_matrix %>%
  group_by(sample_type, group) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate metastatic frequency for each group (BNN and TWT with ecDNA-Positive/Negative)
bnn_neg_freq = as.numeric(round(sample_counts[1,3] / (sample_counts[5,3] + sample_counts[1,3]),2) * 100)
bnn_pos_freq = as.numeric(round(sample_counts[2,3] / (sample_counts[6,3] + sample_counts[2,3]),2) * 100)
twt_neg_freq = as.numeric(round(sample_counts[3,3] / (sample_counts[7,3] + sample_counts[3,3]),2) * 100)
twt_pos_freq = as.numeric(round(sample_counts[4,3] / (sample_counts[8,3] + sample_counts[4,3]),2) * 100)

# Create a frequency data frame for metastatic samples
meta_freq = data.frame(sample = c("TWT_ecDNA-Positive", "TWT_ecDNA-Negative", 
                                  "BNN_ecDNA-Positive", "BNN_ecDNA-Negative"), 
                       freq = c(twt_pos_freq, twt_neg_freq, bnn_pos_freq, bnn_neg_freq))

write.csv(meta_freq, paste0(outdir, "metastatic_frequency.csv"))

#############################################################################################################################

############################################### Copy Number and Gene Analysis #################################################

# Read in GCAP CN (copy number) data and filter for circular ecDNA genes
genes_gcap = read.csv(paste0(outdir,"ALL_gcap_genes_1L.csv"), row.names = 1)
genes_gcap = genes_gcap[genes_gcap$gene_class == "circular", ]  # Keep only circular ecDNA genes

# Create data frames for gene IDs and copy numbers from GCAP
all_genes_gcap = data.frame(sample = genes_gcap$sample, ENSEMBL = genes_gcap$gene_id)
genes_gcap = data.frame(sample = genes_gcap$sample, cn = genes_gcap$total_cn)

# Read in AA CN data and filter for ecDNA classified genes
genes_aa = read.csv(aa_genes)
genes_aa = genes_aa[genes_aa$sample %in% aa_samples & genes_aa$classification == "ecDNA", ]
all_genes_aa = data.frame(sample = genes_aa$sample, SYMBOL = genes_aa$gene)
genes_aa = data.frame(sample = genes_aa$sample, cn = genes_aa$gene_cn)

# Combine GCAP and AA data, then merge with the model matrix to get ecDNA information
cn_data = rbind(genes_aa, genes_gcap)
cn_merge = merge(cn_data, model_matrix, by = "sample")
cn_merge$log2_cn = log2(cn_merge$cn)  # Log-transform copy number values for downstream analysis

# Convert GCAP ENSEMBL IDs to gene symbols
ensembl_ids = all_genes_gcap$ENSEMBL
ensembl_ids = gsub("\\..*", "", ensembl_ids)  # Remove version numbers from ENSEMBL IDs
gene_symbols = AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids,
                                     columns = c("SYMBOL"),
                                     keytype = "ENSEMBL")
all_genes_gcap = na.omit(merge(all_genes_gcap, gene_symbols, by = "ENSEMBL"))  # Remove NA values after merging
all_genes_gcap = all_genes_gcap[-1]  # Drop ENSEMBL IDs column

# Get unique genes from GCAP and AA, then filter for genes present in TWT or BNN subtypes
all_genes = unique(rbind(all_genes_gcap, all_genes_aa))
all_genes = merge(all_genes, model_matrix, by = "sample")
twt_genes = all_genes[all_genes$subtype %in% "TWT", ]
bnn_genes = all_genes[all_genes$subtype %in% "BNN", ]

# Calculate the frequency of genes amplified in each subtype and filter for the most frequent
all_gene_counts = as.data.frame(table(all_genes$SYMBOL))
all_gene_counts = all_gene_counts[order(all_gene_counts$Freq, decreasing = TRUE), ]  # Sort by frequency
write.csv(all_gene_counts, paste0(outdir, "all_ecdna_gene_frequency.csv"))

# Repeat for TWT and BNN
twt_gene_counts = as.data.frame(table(twt_genes$SYMBOL))
twt_gene_counts = twt_gene_counts[order(twt_gene_counts$Freq, decreasing = TRUE), ]  # Sort by frequency
write.csv(twt_gene_counts, paste0(outdir, "twt_ecdna_gene_frequency.csv"))

bnn_gene_counts = as.data.frame(table(bnn_genes$SYMBOL))
bnn_gene_counts = bnn_gene_counts[order(bnn_gene_counts$Freq, decreasing = TRUE), ]  # Sort by frequency
write.csv(bnn_gene_counts, paste0(outdir, "bnn_ecdna_gene_frequency.csv"))

#############################################################################################################################

#################################################### Load Expression Data #####################################################

# Query for RNA-Seq STAR count data from the TCGA-SKCM project
# This pulls gene expression data from the STAR alignment workflow
query = GDCquery(
  project = "TCGA-SKCM",  # Skin Cutaneous Melanoma project in TCGA
  data.category = "Transcriptome Profiling",  # Data category: gene expression
  data.type = "Gene Expression Quantification",  # Quantified gene expression data
  workflow.type = "STAR - Counts"  # STAR alignment tool workflow
)

# Download and prepare the RNA-Seq data
GDCdownload(query)
data = GDCprepare(query)

# Extract the raw counts matrix from the data and adjust sample names
counts = assay(data)
colnames(counts) = substr(colnames(counts), 1, 12)  # Truncate TCGA sample barcodes to 12 characters

# Extract the cell deconvolution TPM counts for further cell type analysis
cell_counts = assay(data, "tpm_unstrand")
colnames(cell_counts) = substr(colnames(counts), 1, 12)  # Ensure matching column names

# Filter the counts to include only samples present in the model matrix
model_matrix = model_matrix[model_matrix$sample %in% colnames(counts), ]
colnames(counts) = make.unique(colnames(counts))  # Ensure unique column names in the counts matrix
counts_filtered = data.frame(counts[, colnames(counts) %in% model_matrix$sample], check.names = FALSE)
cell_counts = data.frame(cell_counts[, colnames(cell_counts) %in% model_matrix$sample], check.names = FALSE)

# Order both the model matrix and the counts matrix by sample name for consistency
model_matrix = model_matrix[order(model_matrix$sample), ]
counts_filtered = counts_filtered[, order(colnames(counts_filtered))]
cell_counts = cell_counts[, order(colnames(cell_counts))]

#############################################################################################################################

################################################## Differential Expression ####################################################

# Conduct differential expression analysis based on group and ecDNA status
# Compare ecDNA positive vs ecDNA negative and group-based differences
ecdna = model_matrix$ecdna  # Define ecDNA status
group = model_matrix$group  # Define experimental group

# Perform DESeq2 differential expression analysis for both ecDNA and group comparisons
dds_ecdna = deseq(counts_filtered, model_matrix, ecdna)
dds_group = deseq(counts_filtered, model_matrix, group)

# Define contrasts for differential expression analysis
# TWT_Positive vs TWT_Negative, BNN_Positive vs BNN_Negative, and ecDNA Positive vs Negative
twt = data.frame(results(dds_group, contrast = c("group", "TWT_ecDNA-Positive", "TWT_ecDNA-Negative")))
bnn = data.frame(results(dds_group, contrast = c("group", "BNN_ecDNA-Positive", "BNN_ecDNA-Negative")))
all = data.frame(results(dds_ecdna, contrast = c("group", "Positive", "Negative")))

# Write the DESeq2 results for each contrast to CSV files for further analysis
write.csv(twt, paste0(outdir, "diff_exp_twt_ecdnaPOS.vs.ecdnaNEG.csv"))
write.csv(bnn, paste0(outdir, "diff_exp_bnn_ecdnaPOS.vs.ecdnaNEG.csv"))
write.csv(all, paste0(outdir, "diff_exp_all_ecdnaPOS.vs.ecdnaNEG.csv"))

#############################################################################################################################

########################################################### GSEA ##############################################################

# Extract significant differentially expressed genes from DESeq2 results (padj ≤ 0.1)
sig_twt = subset(twt, padj <= 0.1)  # Significant genes for TWT contrast
sig_bnn = subset(bnn, padj <= 0.1)  # Significant genes for BNN contrast
sig_all = subset(all, padj <= 0.1)  # Significant genes for all samples (ecDNA Positive vs Negative)

# Conduct Gene Set Enrichment Analysis (GSEA) on the significant DE genes for each contrast
all_gsea = list(gsea(sig_all, "H"), gsea(sig_all, "C1"), gsea(sig_all, "C5"), gsea(sig_all, "C6"),
                gsea(sig_all, "C7"))
# Repeat GSEA for group contrasts
twt_gsea = list(gsea(sig_twt, "H"), gsea(sig_twt, "C1"), gsea(sig_twt, "C5"), gsea(sig_twt, "C6"),
                gsea(sig_twt, "C7"))
bnn_gsea = list(gsea(sig_bnn, "H"), gsea(sig_bnn, "C1"), gsea(sig_bnn, "C5"), gsea(sig_bnn, "C6"),
                gsea(sig_bnn, "C7"))

# Combine the GSEA results into final data frames for each comparison
all_gsea = bind_gsea(all_gsea)  # Combine GSEA results for all samples
twt_gsea = bind_gsea(twt_gsea)  # Combine GSEA results for TWT samples
bnn_gsea = bind_gsea(bnn_gsea)  # Combine GSEA results for BNN samples

# Write the final GSEA results to CSV files
write.csv(all_gsea, paste0(outdir, "all_differential_expression_gsea.csv"))
write.csv(twt_gsea, paste0(outdir, "twt_differential_expression_gsea.csv"))
write.csv(bnn_gsea, paste0(outdir, "bnn_differential_expression_gsea.csv"))

##################################################################################################################################

################################################## Prepare Mixture Matrices ########################################################

# Convert ENSEMBL IDs in the cell count matrix to gene symbols for better interpretation
ensembl_ids = row.names(cell_counts)
ensembl_ids = gsub("\\..*", "", ensembl_ids)  # Remove version numbers from ENSEMBL IDs

# Map ENSEMBL IDs to gene symbols using AnnotationDbi
gene_symbols = AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids,
                                     columns = c("SYMBOL"),
                                     keytype = "ENSEMBL")

# Merge the gene symbols back into the cell count matrix and remove unneeded columns
cell_counts$ENSEMBL = row.names(cell_counts)
cell_counts$ENSEMBL = gsub("\\..*", "", cell_counts$ENSEMBL)
merged_counts = na.omit(merge(gene_symbols, cell_counts, by = "ENSEMBL"))
merged_counts = merged_counts[-c(1, 2)]  # Remove ENSEMBL IDs and SYMBOL columns

# Aggregate counts by gene symbol, taking the mean for duplicated genes
mixture_matrix = merged_counts
mixture_matrix = aggregate(data = mixture_matrix, . ~ Gene, FUN = mean)
row.names(mixture_matrix) = mixture_matrix$Gene  # Set row names to gene symbols
mixture_matrix = mixture_matrix[-1]  # Remove the gene name column

#############################################################################################################################

################################################### Cell Deconvolution ########################################################
# Calculate Tumor Mutational Burden (TMB) for each sample group using somatic mutation data
query = GDCquery(project = "TCGA-SKCM",
                 data.category = "Simple Nucleotide Variation",  # SNP and small indels
                 data.type = "Masked Somatic Mutation",
                 workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query)  # Download the mutation data
maf_file = GDCprepare(query)  # Prepare the mutation data into a MAF (Mutation Annotation Format) object

# Calculate TMB (mutations per megabase) from the MAF file
maf_object = read.maf(maf = maf_file)
tmb_results = tmb(maf = maf_object)
tmb_res = data.frame(sample = substr(tmb_results$Tumor_Sample_Barcode, 1, 12),  # Truncate sample barcodes
                     tmb = tmb_results$total_perMB,  # TMB in mutations per megabase
                     log_tmb = tmb_results$total_perMB_log)  # Log-transformed TMB
tmb_final = merge(tmb_res, model_matrix, by = "sample")  # Merge TMB results with model matrix

# Perform cell deconvolution using xCell to estimate immune and stromal cell types
xcell = deconv(mixture_matrix, "xcell", xcell)  # Use the xCell method for deconvolution
xcell_list = filter_subtypes("xcell", xcell)  # Filter subtypes in the xCell output

# Conduct a Wilcox test to compare cell types between different groups
xcell_all = mannU(xcell_list, "all_pos", "all_neg", "xcell")  # Compare all ecDNA positive vs negative
xcell_all$group = "all"
xcell_bnn = mannU(xcell_list, "bnn_pos", "bnn_neg", "xcell")  # Compare BNN ecDNA positive vs negative
xcell_bnn$group = "bnn"
xcell_twt = mannU(xcell_list, "twt_pos", "twt_neg", "xcell")  # Compare TWT ecDNA positive vs negative
xcell_twt$group = "twt"

# Combine the Wilcox test results and extract significant results (padj ≤ 0.1)
all_decon = rbind(xcell_all, xcell_bnn, xcell_twt)
write.csv(all_decon, paste0(outdir, "xcell_res.csv"))
all_sig = all_decon[all_decon$padj <= 0.1, ]

#############################################################################################################################

################################################### Survival Analysis #########################################################

# Read in survival data and merge with model matrix
survival = read.csv(survival_data)  # Load survival data
surv_merge = merge(survival, model_matrix, by = "sample")  # Merge survival and model data

# Fit Kaplan-Meier survival curves for each group using the survival data
surv_object = Surv(time = survival$time, event = survival$censored)  # Define survival object
km_fit = survfit(surv_object ~ group, data = surv_merge)  # Fit survival model by group

#############################################################################################################################

