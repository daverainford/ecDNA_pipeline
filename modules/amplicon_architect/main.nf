nextflow.enable.dsl=2
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
params.outdir = "/scratch/drainford/skcm_ecdna/ecDNA/results"
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
process AmpliconSuite {
    tag "${sample_id}"
    publishDir "${params.outdir}/", mode: 'copy'

    input:
    tuple val(sample_id), path(mod_bams)

    output:
    path("*.tar.gz"), emit: aa_output

    script:
    def (mod_bam1, mod_bam2) = mod_bams

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
        --bam ${mod_bam2} \\
        --normal_bam ${mod_bam1} \\
        --run_AA \\
        --run_AC 

    find ./ -type f ! -name "*.tar.gz" -exec rm -f {} +
    """
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//Deprecated
// process DataPrep {
//     container '/scratch/drainford/skcm_ecdna/containers/samtools.sif'
    
//     tag "${sample_id}"
//     publishDir "${params.mod_data}", mode: 'move'

//     input:
//     tuple val(sample_id), path(bams), path(bed)

//     output:
//     path "*modified*", emit: modifiedBams

//     script:
//     def (bam1, bam2) = bams

//     """
//     {
//     bam_files=(${bams.join(' ')})
    
//     for file in "\${bam_files[@]}"
//     do
//         outfile=${sample_id}_modified_\${file:5:1}.bam
//         samtools view -@ 1 -H \$file | grep -Fwf <(cut -f1 ${bed} | sort -u) - | \
//         samtools reheader - <(samtools view -L ${bed} -b \$file) > \$outfile &
//     done
//     wait

//     for bam in *.bam
//     do
//         samtools index -@ 1 -b \${bam} &
//     done
//     wait
//     } > ./${sample_id}.log 2>&1
//     """
// }

///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
