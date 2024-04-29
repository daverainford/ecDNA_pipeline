nextflow.enable.dsl=2

// Define parameters
params.data = "./data/test_fastq/" 
params.outdir = "./results/test_results/"
params.user = "HPC"
params.test = false
params.help = false

// Define process to test the AmpliconSuite pipeline 
process TestInstall {
    tag "${sample_id}"
    publishDir "${params.outdir}/${sample_id}/", mode: 'move'

    input:
    tuple val(sample_id), path(reads)

    output:
    path("*.tar.gz")

    script:
    // Split FASTQ files into read pairs
    def (read1, read2) = reads

    // Pull container, run AmpliconSuite test, clean up work dir.
    """
    module load singularity
    module load git

    export HOME=/home/${params.user}/

    singularity pull ./ampliconsuite-pipeline.sif library://jluebeck/ampliconsuite-pipeline/ampliconsuite-pipeline:1.3.1
    git clone https://github.com/AmpliconSuite/AmpliconSuite-pipeline

    ./AmpliconSuite-pipeline/singularity/run_paa_singularity.py \\
        -s ${sample_id} \\
        -t 5 \\
        -o ./ \\
        --fastqs ${read1} ${read2} \\
        --ref GRCh38 \\
        --run_AA \\
        --run_AC &> ./testInstall.log

    find ./ -type f ! -name "*.tar.gz" -exec rm -f {} +
    """
}
