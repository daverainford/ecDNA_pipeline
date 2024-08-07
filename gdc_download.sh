#!/bin/bash
#SBATCH --job-name=gdc_client
#SBATCH --mail-type=ALL     
#SBATCH --mail-user=drainford@tgen.org  
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=40                                    
#SBATCH --time=96:00:00               
#SBATCH --output=/scratch/drainford/logs/gdc_client_rna.log

# Load singularity
module load singularity

# Download TCGA data to data download directory
singularity exec /scratch/drainford/containers/gdc-client_1.6.1.sif gdc-client download \
    -n 40 \
    -m download_data/rna_manifest.txt \
    -t download_data/*token* \
    -d /scratch/drainford/ecDNA_pipeline/data/rna/

# Empty all subdirectories in data download directory
find /scratch/drainford/ecDNA_pipeline/data/ -mindepth 2 -type f -exec mv -t /scratch/drainford/ecDNA_pipeline/data/ {} +

# Delete the now empty directories
find /scratch/drainford/ecDNA_pipeline/data/ -mindepth 1 -type d -exec rm -r {} +

# Remove uneeded index and parcel files
rm /scratch/drainford/ecDNA_pipeline/data/*.bai /scratch/drainford/ecDNA_pipeline/data/*.parcel