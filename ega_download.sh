#!/bin/bash
#SBATCH --job-name=ega_client
#SBATCH --mail-type=ALL     
#SBATCH --mail-user=drainford@tgen.org  
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=40                                    
#SBATCH --time=96:00:00               
#SBATCH --output=/scratch/drainford/logs/ega_client.log

module load singularity

# Store EGA file IDs in variable from batch file
files=$(cat /scratch/drainford/ecDNA_pipeline/download_data/pcawg_batch1.txt)

# Iterate over file ids to download BAMs
for file in ${files}
do
    singularity exec /scratch/drainford/containers/pyega3_5.2.0--pyhdfd78af_0.sif pyega3 \
        -c 40 \
        -cf /scratch/drainford/ecDNA_pipeline/download_data/credentials.json \
        fetch ${file:0:15} \
        --output-dir /scratch/drainford/ecDNA_pipeline/mela_au_wgs_bams/ 
done
