#!/bin/bash

#SBATCH --job-name=vpt1-array
#SBATCH -A {account}
#SBATCH -p standard
#SBATCH -n 6
#SBATCH --array=0-11
#SBATCH -t 06:00:00
#SBATCH --mem=60GB
#SBATCH --output=/scratch/{user}/2024-01-23_MERFISH/slurm-scripts/slurm-out/vpt1-%a.out

# Setup
dirs=("E003-VPT017" "E003-VPT018" "E003-VPT019" "E003-VPT020" "E003-VPT021" "E003-VPT022" "E003-VPT023" "E003-VPT024" "E003-VPT025" "E003-VPT026" "E003-VPT027" "E003-VPT028")

dir=${dirs[$SLURM_ARRAY_TASK_ID]}

module purge
module load anaconda
source activate vpt_env

echo "Processing $dir..."

vpt --verbose --processes 4 run-segmentation \
--segmentation-algorithm "/scratch/{user}/2024-01-23_MERFISH/segmentation-algorithms/${dir}.json" \
--input-images="/scratch/{user}/2024-01-23_MERFISH/E003-original/region_0/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
--input-micron-to-mosaic /scratch/{user}/2024-01-23_MERFISH/E003-original/region_0/images/micron_to_mosaic_pixel_transform.csv \
--output-path "/scratch/{user}/2024-01-23_MERFISH/${dir}/analysis_outputs" \
--tile-size 2400 \
--tile-overlap 200

# Capture the exit status of the last command
status=$?

# Check if the vpt command was successful
if [ $status -eq 0 ]; then
    echo "Processing of $dir completed successfully."
else
    echo "Processing of $dir failed with status code $status."
fi
