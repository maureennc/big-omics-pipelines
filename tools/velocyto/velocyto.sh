#!/bin/bash

#SBATCH --job-name=naive-velocyto
#SBATCH -A harrislab
#SBATCH -p standard
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH -t 06:00:00
#SBATCH --mem=256GB

module purge
module load anaconda
module load samtools

samtools --version

source activate velocyto

echo "Processing naive sample..."

velocyto run10x --samtools-threads 8 \
                --samtools-memory 120000 \
                "/scratch/mnc3ra/cite-seq/velocyto/naive_cellranger_output_copy" \
                /scratch/mnc3ra/cite-seq/velocyto/genes.gtf