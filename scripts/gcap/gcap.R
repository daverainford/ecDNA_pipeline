# Add arguments to execute gcap.R from the command line
###########################################################Argparse###########################################################
library(argparse)

# Add and define the --ascat_files and --outdir arguments
parser <- ArgumentParser()

parser$add_argument('--ascat_files', required=TRUE, help="Path to directory containing ASCAT CNV segments TXT files.")
parser$add_argument('--outdir', required=TRUE, help="Path to output directory.")
args = parser$parse_args()
ascat_cnv = args$ascat_files
outdir = args$outdir

if (!dir.exists(ascat_cnv)) {
    stop(paste("This directory does not exist:", ascat_cnv))
}

# Check if gdc_data directory exists
if (!dir.exists(gdc_data)) {
  stop(paste("This directory does not exist:", gdc_data))
}
###############################################################################################################################

# Run GCAP workflow to identify ecDNA from Affy genotyping chip CNV segments
############################################################ GCAP #############################################################
library(gcap)

# List all files in the ascat directory for processing
ascat_files = list.files(ascat_cnv)  # Get the list of ASCAT output files

# Read the first ASCAT file and create a dataframe with relevant columns
ascat = read.table(paste0(ascat_cnv, ascat_files[1]), sep = "\t", header = TRUE)
ascat = data.frame(chromosome = ascat$Chromosome,  # Chromosome number
                   start = ascat$Start,  # Start position of the segment
                   end = ascat$End,  # End position of the segment
                   total_cn = ascat$Copy_Number,  # Total copy number of the segment
                   minor_cn = ascat$Minor_Copy_Number,  # Minor allele copy number
                   sample = substr(ascat_files[1], 1, 12))  # Extract sample ID from file name

# Loop through the rest of the ASCAT files and append their data to the existing dataframe
for (file in ascat_files[2:length(ascat_files)]) {
  data = read.table(paste0(ascat_cnv, file), sep = "\t", header = TRUE)  # Read each ASCAT file
  df = data.frame(chromosome = data$Chromosome,  # Create dataframe for the current file
                  start = data$Start,
                  end = data$End,
                  total_cn = data$Copy_Number,
                  minor_cn = data$Minor_Copy_Number,
                  sample = substr(file, 1, 12))  # Extract sample ID from file name
  ascat = rbind(ascat, df)  # Append the data to the main ASCAT dataframe
}

# Write the combined ASCAT data to a CSV file
write.csv(ascat, paste0(outdir, "gcap_data.csv"))  # Save the full ASCAT data in CSV format

# Perform GCAP ASC-N workflow on the combined ASCAT data
# "XGB11" is the model being used, and the genome build is hg38
gcap_res = gcap.ASCNworkflow(ascat, 
                             outdir = tempdir(),  # Use a temporary directory for output
                             model = "XGB11",  # Use the XGB11 model for analysis
                             genome_build = "hg38",  # Human genome build hg38
                             tightness = 1L)  # Optimal tightness parameter for GCAP analysis

# Extract the sample summary and gene summary from the GCAP results
res_summary = gcap_res$sample_summary  # Get the summary of samples analyzed
gene_summary = gcap_res$data  # Get the gene-level summary from the analysis

# Write the sample summary and gene summary to CSV files for further analysis
write.csv(res_summary, paste0(outdir,"ALL_gcap_res_1L.csv"))  # Save the sample summary
write.csv(gene_summary, paste0(outdir, "ALL_gcap_genes_1L.csv"))  # Save the gene-level summary
#############################################################################################################################
