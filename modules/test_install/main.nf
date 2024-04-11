nextflow.enable.dsl=2
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
params.outdir = "/scratch/drainford/skcm_ecdna/ecDNA/results"
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
process TestInstall {
    tag "${sample_id}"
    publishDir "${params.outdir}/test/", mode: 'move'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.tar.gz"), emit: aa_output

    script:
    def (read1, read2) = reads

    """
    module load singularity
    module load git

    export AA_DATA_REPO=/scratch/drainford/skcm_ecdna/aa_data_repo/
    export HOME=/home/drainford/

    singularity pull ./ampliconsuite-pipeline.sif library://jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline:1.3.1
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline

    mkdir ./logs/
    ./AmpliconSuite-pipeline/singularity/run_paa_singularity.py \\
        -s ${sample_id} \\
        -t 2 \\
        -o ${params.outdir}/test/ \\
        --fastqs ${read1} ${read2} \\
        --ref hg19 \\
        --run_AA \\
        --run_AC &> ./testInstall.log

    find ./ -type f ! -name "*.tar.gz" -exec rm -f {} +
    """
}
