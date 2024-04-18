nextflow.enable.dsl=2

// Define params. This will be removed, as the params are not currently carrying over from main.nf. Need to debug.
params.outdir = "/scratch/drainford/skcm_ecdna/ecDNA/results"

// Define process to execute ampliconSuite pipeline.
process AmpliconSuite {
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

    export AA_DATA_REPO=/home/drainford/aa_data_repo/
    export HOME=/home/drainford/

    singularity pull ./ampliconsuite-pipeline.sif library://jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline:1.3.1
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline

    ./AmpliconSuite-pipeline/singularity/run_paa_singularity.py \\
        -s ${sample_id} \\
        -t 2 \\
        -o ./ \\
        --ref GRCh38 \\
        --bam ${normal} \\
        --normal_bam ${tumor} \\
        --run_AA \\
        --run_AC 

    find ./ -type f ! -name "*.tar.gz" -exec rm -f {} +
    """
}
