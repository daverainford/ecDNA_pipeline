# Overview
This is a Nextflow pipeline to handle containerized [Link AmpliconArchitect]([URL](https://github.com/virajbdeshpande/AmpliconArchitect)) inputs and outputs.

# Data Format
This pipeline takes a matched tumor/normal pair of BAM file as input.

# Usage
```
nextflow run main.nf \
  --data <path to BAM files> \
  --outdir <path to output directory> \
  --test <true if you want to test the pipeline; default false>
```
