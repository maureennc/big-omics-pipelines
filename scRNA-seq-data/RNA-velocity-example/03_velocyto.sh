#!/bin/bash

#SBATCH --job-name=INF1-velocyto
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

echo "Processing INF1..."

cd /project/harrislab/mnc3ra/scrna-seq/01-2024_cite-seq/velocyto/

velocyto run10x \
--samtools-threads 8 \
--samtools-memory 120000 \
-m input-files/mm10_rmsk-mito-edit.gtf \
input-files/INF1/INF1_cellranger_output_velocyto \
input-files/mito-edit.genes.gtf
