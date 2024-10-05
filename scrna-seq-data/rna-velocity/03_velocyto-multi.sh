#!/bin/bash

#SBATCH --job-name=velocyto
#SBATCH -A {account}
#SBATCH -p standard
#SBATCH -n 6
#SBATCH --array=0-2
#SBATCH -t 06:00:00
#SBATCH --mem=60GB
#SBATCH --output="/scratch/$USER/scRNA-seq/velocyto/slurm-scripts/slurm-out/velocyto-%a.out"

# Setup
dirs=("naive" "INF1" "INF2")

# Get the directory based on the SLURM_ARRAY_TASK_ID
dir=${dirs[$SLURM_ARRAY_TASK_ID]}

module purge
module load anaconda
source activate scvelo

echo "Processing $dir..."


velocyto run10x "/scratch/$USER/scRNA-seq/velocyto/${dir}_cellranger_output" /scratch/$USER/scRNA-seq/velocyto/genes.gtf
