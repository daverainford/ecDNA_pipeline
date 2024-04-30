//Deprecated process

// nextflow.enable.dsl=2

// // Define parameters
// params.data = "./data/test_fastq/" 
// params.outdir = "./results/test_results/"
// params.user = "HPC"
// params.test = false
// params.help = false

// // Define params. This will be removed, as the params are not currently carrying over from main.nf. Need to debug.
// params.outdir = "/scratch/drainford/skcm_ecdna/ecDNA/results"

// // Define process to sumarrize AmpliconSuite results
// process DataSummary {
//     publishDir "${params.outdir}/", mode: 'move'

//     input:
//     path(aa_output)

//     output:
//     path("*table.tsv"), emit: summary

//     // Extract *results_table.tsv from AmpliconSuite output *.tar.gz
//     script:
//     """
//     {
//     tar -zxf *.tar.gz 
//     mv *classification/*table.tsv .
//     find ./ -type f ! -name "*table.tsv" -exec rm -f {} +
//     } >> DataSummary.log 2>&1
//     """
// }
// ///////////////////////////////////////////////////////////////////////////////////////////////////
// ///////////////////////////////////////////////////////////////////////////////////////////////////