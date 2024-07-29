nextflow.enable.dsl=2

// Load modules
include {TestInstall} from './modules/test_install/main.nf'
include {AmpliconSuite} from './modules/amplicon_architect/main.nf'
include {Gcap} from './modules/gcap/main.nf'

// Define parameters
params.data = "./data/test_fastq/" 
params.outdir = "./results/test_results/"
params.user = "HPC"
params.gcap_conatiner = "/scratch/drainford/containers/gcap_latest.sif"
params.wes = false
params.wgs = false
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
--user      HPC user account ID
--test      Argument to run TestInstall process.

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
    } 

    if (params.wgs) {
        def bams_channel = Channel.fromFilePairs("${params.data}/*_{T,N}.bam")
            .ifEmpty {throw new RuntimeException("No BAM files found matching pattern in ${params.data}/")}
        def ampliconOutputs = AmpliconSuite(bams_channel)
    }

    if (params.wes) {
        def bams_channel = Channel.fromFilePairs("${params.data}/*_{T,N}.bam")
            .ifEmpty {throw new RuntimeException("No BAM files found matching pattern in ${params.data}/")}
        def gcapOutputs = Gcap(bams_channel)
    }
}
