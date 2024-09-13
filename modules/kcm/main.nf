nextflow.enable.dsl=2

// Define parameters
params.data = "./data/test_fastq/" 
params.outdir = "./results/test_results/"
params.user = "HPC"
params.test = false
params.help = false

// Define process to test the AmpliconSuite pipeline 
process Kmc {
    
}