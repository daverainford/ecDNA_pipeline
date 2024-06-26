#!/bin/bash
#SBATCH --job-name=gdc_client
#SBATCH --mail-type=ALL     
#SBATCH --mail-user=drainford@tgen.org  
#SBATCH --ntasks=1 
#SBATCH --cpus-per-task=20                                    
#SBATCH --time=36:00:00               
#SBATCH --output=/scratch/drainford/logs/gdc_client.log

# Load singularity
module load singularity

# Download TCGA data to data download directory
singularity exec /scratch/drainford/containers/gdc-client_1.6.1.sif gdc-client download \
    -n 20 \
    -m download_data/*batch* \
    -t download_data/*token* \
    -d /scratch/drainford/ecDNA/data/

# Empty all subdirectories in data download directory
find /scratch/drainford/ecDNA/data/ -mindepth 2 -type f -exec mv -t /scratch/drainford/ecDNA/data/ {} +

# Delete the now empty directories
find /scratch/drainford/ecDNA/data/ -mindepth 1 -type d -exec rm -r {} +

# Remove uneeded index and parcel files
rm /scratch/drainford/ecDNA/data/*.bai /scratch/drainford/ecDNA/data/*.parcel

# # Get md5sums for bams all files in download directory
# md5sum /scratch/drainford/ecDNA/data/* > checksums.txt