#!/bin/bash

#SBATCH --job-name=vpt2-array
#SBATCH -A {account}
#SBATCH -p standard
#SBATCH -n 6
#SBATCH --array=0-11
#SBATCH -t 03:00:00 
#SBATCH --mem=60GB
#SBATCH --output=/scratch/{user}/2024-01-23_MERFISH/slurm-scripts/slurm-out/vpt2-%a.out

# Setup
dirs=("E003-VPT017" "E003-VPT018" "E003-VPT019" "E003-VPT020" "E003-VPT021" "E003-VPT022" "E003-VPT023" "E003-VPT024" "E003-VPT025" "E003-VPT026" "E003-VPT027" "E003-VPT028")

dir=${dirs[$SLURM_ARRAY_TASK_ID]}

module purge
module load anaconda
source activate vpt_env

echo "Processing $dir..."

vpt --verbose partition-transcripts \
--input-boundaries "/scratch/{user}/2024-01-23_MERFISH/${dir}/analysis_outputs/cellpose_micron_space.parquet" \
--input-transcripts "/scratch/{user}/2024-01-23_MERFISH/E003-original/region_0/detected_transcripts.csv" \
--output-entity-by-gene "/scratch/{user}/2024-01-23_MERFISH/${dir}/analysis_outputs/${dir}_cell_by_gene.csv" \
--output-transcripts "/scratch/{user}/2024-01-23_MERFISH/${dir}/analysis_outputs/${dir}_detected_transcripts.csv"

# Capture the exit status of the last command
status=$?

# Check if the vpt command was successful
if [ $status -eq 0 ]; then
    echo "Partitioning of transcripts for $dir completed successfully."
else
    echo "Partitioning of transcripts for $dir failed with status code $status."
fi
