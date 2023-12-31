#!/bin/bash

#SBATCH -n 6
#SBATCH -t 0:30:00
#SBATCH -p standard
#SBATCH -A harrislab
#SBATCH --mem=20GB

module purge
module load fastqc

fastqc /project/harrislab/naive_vs_infected_whole_brain/trimmed_paired/I2_R1_001_trimmed.fastq.gz -o /project/harrislab/naive_vs_infected_whole_brain/qc

fastqc /project/harrislab/naive_vs_infected_whole_brain/trimmed_paired/I2_R2_001_trimmed.fastq.gz -o /project/harrislab/naive_vs_infected_whole_brain/qc
