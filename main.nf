nextflow.enable.dsl=2

// Load modules
include {TestInstall} from './modules/test_install/main.nf'
include {AmpliconSuite} from './modules/amplicon_architect/main.nf'
include {DataSummary} from './modules/data_summary/main.nf'

// Define parameters
params.data = "./data/test_fastq/" 
params.outdir = "./results/test_results/"
params.user = "HPC"
params.test = false
params.help = false

// Help message
def helpMessage = 
"""
Usage:
nextflow run script.nf --data <path> --outdir <path> --test <Boolean>

Options:
--data      Path to the directory containing input BAM files OR test FASTQ files.
--outdir    Path to the directory where AmpliconArchitect with publish the results.
--test      Boolean to run the TestInstall process. Default is false.

--help      Print this help message.
"""

// If "--help true" is specified, print help message
if (params.help) {
    println(helpMessage)
    System.exit(0) 
}

// Define nextflow workflow and data channels
workflow {
    // If "--test true" is specified, test pipeline, else run as normal
    if (params.test) {
        def test_reads_channel = Channel.fromFilePairs("${params.data}/*_{1,2}.fastq")
            .ifEmpty {throw new RuntimeException("No FASTQ files found matching pattern in ${params.data}/")}
        TestInstall(test_reads_channel)
    } else {
        def bams_channel = Channel.fromFilePairs("${params.data}/*_{T,N}.bam")
            .ifEmpty {throw new RuntimeException("No BAM files found matching pattern in ${params.data}/")}
        def ampliconOutputs = AmpliconSuite(bams_channel)
        DataSummary(ampliconOutputs.aa_output)
    }
}
