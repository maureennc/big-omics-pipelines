#!/bin/bash

#SBATCH --job-name=INF1-samtools
#SBATCH -A harrislab
#SBATCH -p standard
#SBATCH -n 1
#SBATCH --cpus-per-task=8
#SBATCH -t 09:00:00
#SBATCH --mem=100GB

module purge
module load samtools

samtools --version

samtools sort -t CB -O BAM -@ 7 \
	-o /project/harrislab/mnc3ra/scrna-seq/01-2024_cite-seq/samtools/output-files/INF1_cellsorted_genome_bam.bam \
	/project/harrislab/mnc3ra/scrna-seq/01-2024_cite-seq/samtools/input-files/INF1_possorted_genome_bam.bam
