# Overview
This is a Nextflow pipeline to handle containerized [AmpliconSuite](https://github.com/AmpliconSuite/AmpliconSuite-pipeline) inputs and outputs.

# Data Format
This pipeline takes a matched tumor/normal pair of BAM file as input. 
It is reccomended to test AmpliconSuite upon initial install, which can be done here by: ```--test``` flag.
The AmpliconArchitect test FASTQ files can be found by following the testing instructions in the link above.

# Usage
```
nextflow run main.nf \
  --user HPC user account ID
  --data <path to BAM or test FASTQ files> \
  --outdir <path to output directory> \
  --test <pass if you want to run the TestInstall process; default false>
  -bg Use this flag to run a batch job on SLURM
```
