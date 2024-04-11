nextflow.enable.dsl=2
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
params.outdir = "/scratch/drainford/skcm_ecdna/ecDNA/results"
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
process DataSummary {
    publishDir "${params.outdir}/", mode: 'move'

    input:
    path(aa_output)

    output:
    path("*table.tsv"), emit: summary

    script:
    """
    {
    tar -zxf *.tar.gz 
    mv *classification/*table.tsv .
    find ./ -type f ! -name "*table.tsv" -exec rm -f {} +
    } >> DataSummary.log 2>&1
    """
}
///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////