#!/bin/bash

#SBATCH -n 6
#SBATCH -t 0:30:00
#SBATCH -p standard
#SBATCH -A harrislab
#SBATCH --mem=20GB

module purge

module load trimmomatic

java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads 5 -phred33 /project/harrislab/naive_vs_infected_whole_brain/fastq-files/I2_R1_001.fastq.gz /project/harrislab/naive_vs_infected_whole_brain/fastq-files/I2_R2_001.fastq.gz /project/harrislab/naive_vs_infected_whole_brain/trimmed_paired/I2_R1_001_trimmed.fastq.gz /project/harrislab/naive_vs_infected_whole_brain/trimmed_unpaired/I2_R1_001_unpaired.fastq.gz /project/harrislab/naive_vs_infected_whole_brain/trimmed_paired/I2_R2_001_trimmed.fastq.gz /project/harrislab/naive_vs_infected_whole_brain/trimmed_unpaired/I2_R2_001_unpaired.fastq.gz ILLUMINACLIP:TruSeq2-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36
