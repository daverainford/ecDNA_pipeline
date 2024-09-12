# Branch Overview 
The all branch of ecDNA_pipeline contains all pertinent code for any analyses related to the SKCM-ecDNA project.
Below is a step by step procedure on how to run all code and reproduce analyses 

# Data Download
The below workflows take either WGS Bam files or ASCAT3 allele specific copy number segments.
These can be download from the TCGA/EGA with the included helper script gdc_download.sh and ega_download.sh
You will have to modify the scripts to download the data you are specifically after.

# Amplicon Architect Overview
ecDNA_pipeline/main.nf is a Nextflow pipeline to handle containerized [AmpliconSuite](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) inputs and outputs.

# Amplicon Architect Data Format
This pipeline takes a matched tumor/normal pair of BAM file as input. 
It is reccomended to test AmpliconSuite upon initial install, which can be done here by: ```--test``` flag.
The AmpliconArchitect test FASTQ files can be found by following the testing instructions in the link above.

# Amplicon Architect Usage
```
nextflow run main.nf \
  --user HPC user account ID
  --data <path to BAM or test FASTQ files> \
  --outdir <path to output directory> \
  --test <true if you want to test the pipeline; default false>
  -bg Use this flag to run a batch job on SLURM
```

# GCAP Overview
ecDNA_pipeline/R/gcap.R is a script to run GCAP [GCAP](https://github.com/ShixiangWang/gcap).

# GCAP Data Format
This script takes allele-specific copy number segments as input.
It outputs two CSV files: one with classified samples and one with gene level amplification/copy-number information
The script is currently only structured to work with ASCAT3 segment outputs.

# GCAP Usage
```
Rscript gcap.R \
  --ascat_files <path to directory containing ASCAT3 copy number segments> \
  --outdir <path to the directory to write outputs to>
```

# main.R Overview
ecDNA_pipeline/R/gcap.R is an R script that runs several analyses downstream of ecDNA classifiction including:
ecDNA frequency calculation
sample type frequency calculation
ecDNA gene-level amplification frequency and copy number calculations
DESeq2 differential expression
GSEA
XCell cell deconvolution
survival analysis.

# main.R Data Format
main.R takes three CSV data files as input: a model matrix defining group comparisons, a gene summary generated
from Amplicon Architect, and sample-level survival data donwload from the TCGA analysis center.

# main.R Usage
```
Rscript main.R \
  --model <path to model matrix CSV> \
  --aa_genes <path to AmpliconArchitect gene summary CSV> \
  --survival_data <path to survival data CSV>
  --outdir <path to the directory to write outputs to>
```