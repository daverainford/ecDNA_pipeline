# Overview
This is a Nextflow pipeline to handle containerized [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect) inputs and outputs.

# Data Format
This pipeline takes a matched tumor/normal pair of BAM file as input. 
It is reccomended to test AmpliconArchitect upon initial install, which can be done here by: ```--test``` flag.
The AmpliconArchitect test FASTQ files can be found by following the testing instructions in the link above.

# Usage
```
nextflow run main.nf \
  --user HPC user account ID
  --data <path to BAM or test FASTQ files> \
  --outdir <path to output directory> \
  --test <true if you want to test the pipeline; default false>
  -bg Use this flag to run a batch job on SLURM
```
