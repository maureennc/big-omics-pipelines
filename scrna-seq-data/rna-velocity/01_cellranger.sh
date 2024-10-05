#!/bin/bash

#SBATCH --job-name=INF1-mito-count
#SBATCH -A harrislab
#SBATCH -p standard
#SBATCH -c 16
#SBATCH -t 12:00:00
#SBATCH --mem=150GB

module purge
module load cellranger

# Navigate to output directory
cd /project/harrislab/mnc3ra/scrna-seq/01-2024_cite-seq/cellranger/output-files-velocyto

cellranger count --id=INF1_cellranger_output_velocyto \
--libraries=../input-files/libraries/library-INF1.csv \
--transcriptome=../input-files/mm10/mkref-mm10-mito-edit \
--feature-ref=../input-files/feature-ref/feature-ref.csv \
--localcores=16 \
--localmem=128 \
