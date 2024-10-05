#!/bin/bash

#SBATCH -n 7
#SBATCH -t 05:00:00
#SBATCH -p standard
#SBATCH -A {account}
#SBATCH --mem=75G

module purge
module load gcc/9.2.0 salmon/1.5.1

salmon quant -i /project/harrislab/salmon_mm10_index -p 6 -l A -1 /project/harrislab/naive_vs_infected_whole_brain/trimmed_paired/I2_R1_001_trimmed.fastq.gz -2 /project/harrislab/naive_vs_infected_whole_brain/trimmed_paired/I2_R2_001_trimmed.fastq.gz --validateMappings -o I2_salmon
