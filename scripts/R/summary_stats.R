# Define groups for experimental comaprisons
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

# Calculate ecDNA frequency for each group
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

# Calculate percent of metastatic samples for each group
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

# Calculate ecDNA amplification gene frequency and gene-level ecDNA copy number
############################################### Copy Number and Gene Analysis #################################################
# Read in GCAP CN (copy number) data and filter for circular ecDNA genes
genes_gcap = read.csv(gcap_genes, row.names = 1)
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
write.csv(cn_merge, paste0(outdir, "ecdna_gene_level_cn.csv"))


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

# Fit survival analysis Kaplan-Meier model
################################################### Survival Analysis #########################################################
# Read in survival data and merge with model matrix
survival = read.csv(survival_data)  # Load survival data
surv_merge = merge(survival, model_matrix, by = "sample")  # Merge survival and model data

# Fit Kaplan-Meier survival curves for each group using the survival data
surv_object = Surv(time = survival$time, event = survival$censored)  # Define survival object
km_fit = survfit(surv_object ~ group, data = surv_merge)  # Fit survival model by group
#############################################################################################################################