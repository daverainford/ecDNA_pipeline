nextflow.enable.dsl=2

// Define parameters
params.data = "./data/test_fastq/" 
params.outdir = "./results/test_results/"
params.user = "HPC"
params.test = false
params.help = false

// Define process to execute ampliconSuite pipeline.
process Gcap {
    tag "${sample_id}"
    publishDir "${params.outdir}/", mode: 'move'

    input:
    tuple val(sample_id), path(mod_bams)

    output:
    path("*.tar.gz"), emit: aa_output

    script:
    // Seperate T/N bams for use in pipeline.
    def (normal, tumor) = mod_bams

    // Pull container, run AmpliconSuite, clean up work dir.
    """
    module load git
    module load singularity

    export HOME=/home/${params.user}/
    
    singularity pull ./ampliconsuite-pipeline.sif library://jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline:1.3.1
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline

    ./AmpliconSuite-pipeline/singularity/run_paa_singularity.py \\
        -s ${sample_id} \\
        -t 5 \\
        -o ./ \\
        --ref GRCh38 \\
        --bam ${tumor} \\
        --normal_bam ${normal} \\
        --run_AA \\
        --run_AC 

    find ./ -type f ! -name "*.tar.gz" -exec rm -f {} +
    """
}
