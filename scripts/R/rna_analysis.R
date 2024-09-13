# Load count matricies for both differential expression and cell deconvolution
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
GDCdownload(query, directory = gdc_data)
data = GDCprepare(query, directory = gdc_data)

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

# Perform differential expression analyses on all groups
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

# Perform GSEA on each groups differentail expression results
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

# Prepare the mixture matrix (count matrix) for cell deconvolution input
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
merged_counts = merged_counts[-1]  # Remove ENSEMBL IDs and SYMBOL columns

# Aggregate counts by gene symbol, taking the mean for duplicated genes
mixture_matrix = merged_counts
mixture_matrix = aggregate(data = mixture_matrix, . ~ SYMBOL, FUN = mean)
row.names(mixture_matrix) = mixture_matrix$Gene  # Set row names to gene symbols
mixture_matrix = mixture_matrix[-1]  # Remove the gene name column

#############################################################################################################################

# Calculate TMB and perform cell deconvolution with xcell
################################################### Cell Deconvolution ########################################################
# Calculate Tumor Mutational Burden (TMB) for each sample group using somatic mutation data
query = GDCquery(project = "TCGA-SKCM",
                 data.category = "Simple Nucleotide Variation",  # SNP and small indels
                 data.type = "Masked Somatic Mutation",
                 workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(query, directory = gdc_data)  # Download the mutation data
maf_file = GDCprepare(query, directory = gdc_data)  # Prepare the mutation data into a MAF (Mutation Annotation Format) object

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