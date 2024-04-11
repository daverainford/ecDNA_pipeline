# Overview
This is a Nextflow pipeline to handle containerized [AmpliconArchitect](https://github.com/virajbdeshpande/AmpliconArchitect) inputs and outputs.

# Data Format
This pipeline takes a matched tumor/normal pair of BAM file as input. 
It is reccomended to test AmpliconArchitect upon initial install, which can be done here by: ```--test true```.

# Usage
```
nextflow run main.nf \
  --data <path to BAM or test FASTQ files> \
  --outdir <path to output directory> \
  --test <true if you want to test the pipeline; default false>
```
