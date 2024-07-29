nextflow.enable.dsl=2

// Define parameters
params.data = "./data/test_fastq/" 
params.outdir = "./results/test_results/"
params.user = "HPC"
params.gcap_conatiner = "/scratch/drainford/containers/gcap_latest.sif"
params.wes = false
params.wgs = false
params.test = false
params.help = false

// Define process to execute ampliconSuite pipeline.
process Gcap {
    tag "${sample_id}"
    publishDir "${params.outdir}/", mode: 'move'
    container params.gcap_container

    input:
    tuple val(sample_id), path(mod_bams), path("gcap.R")

    output:
    path("*.tsv"), emit: aa_output

    script:
     // Pull container, run GCAP, clean up work dir.
    """
    Rscript gcap.R
    find ./ -type f ! -name "*.tsv" -exec rm -f {} +
    """
}
