#!/bin/bash
#SBATCH --job-name=ega_client
#SBATCH --mail-type=ALL     
#SBATCH --mail-user=drainford@tgen.org  
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20                                    
#SBATCH --time=96:00:00               
#SBATCH --output=/scratch/drainford/logs/ega_client.log

module load singularity

singularity exec pyega3_5.2.0--pyhdfd78af_0.sif pyega3 \
    -c 20 \
    -cf /scratch/drainford/ecDNA_pipeline/download_data/credentials.json \
    fetch EGAD00001002157 \
    --output-dir /scratch/drainford/ecDNA_pipeline/
