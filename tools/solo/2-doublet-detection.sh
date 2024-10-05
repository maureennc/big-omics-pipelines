#!/bin/bash

#SBATCH -A lukenslab
#SBATCH -p standard
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=solo_%j.log

module purge
module load anaconda

source activate sc-seq

python3 /scratch/mnc3ra/tbi_snseq/solo/python/2-doublet-detection.py
