# Add arguments for execution fo main.R from the command line
###########################################################Argparse###########################################################

library(argparse)

# Add command line arguments with argparse
parser <- ArgumentParser()

parser$add_argument('--model', required=TRUE, help="Path to the model matrix CSV. 
                    Required cols: sample, ecdna, subtype, group, sample_type, method.")
parser$add_argument('--aa_genes', required=TRUE, help="Path to the AmpliconArchitect gene summary CSV.")
parser$add_argument('--survival_data', required=TRUE, help="Path to the survival_data CSV.
                    Required cols: sample, time, censored.")
parser$add_argument('--gdc_data', required=TRUE, help="Path to directory containing TCGA transcriptome/MAF data, or directory you want to download it in.
                                                        If you downloading it for the first time, it is important you do it in the same directory as the main.R
                                                        script.")
parser$add_argument('--gcap_genes', required=TRUE, help="Path to ALL_gcap_genes_1L.csv.")
parser$add_argument('--outdir', required=TRUE, help="Path to output directory.")

# Store command line arguments as variables
args = parser$parse_args()
model = args$model
aa_genes = args$aa_genes
survival_data = args$survival_data
gdc_data = normalizePath(as.character(args$gdc_data), mustWork = FALSE)
gcap_genes = args$gcap_genes
outdir = normalizePath(as.character(args$outidr), mustWork = FALSE)

# Source functions for use in analysis
source(file.path(getwd(), "utils.R"))

# Check if files/directories provided in arguments exist
file_check(model)
file_check(aa_genes)
file_check(survival_data)
file_check(gcap_genes)
dir_check(gdc_data)
###############################################################################################################################

# Load all libraries needed for analyses
###########################################################Load Libraries###########################################################

library(GenomicFeatures)
library(org.Hs.eg.db)
library(immunedeconv)
library(tidyverse)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(maftools)
library(msigdbr)
library(survival)
library(survminer)
library(DESeq2)
library(gcap)
library(clusterProfiler)

# Execute summary_stats.R
source(file.path(getwd(), "summary_stats.R"))

# Execute rna_analysis
source(file.path(getwd(), "rna_analysis.R"))