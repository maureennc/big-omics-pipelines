#!/bin/bash

#SBATCH --job-name=vpt1
#SBATCH -A {account}
#SBATCH -p standard
#SBATCH -n 6
#SBATCH -t 06:00:00
#SBATCH --mem=40GB

module purge
module load anaconda

source activate vpt_env

vpt --verbose --processes 4 run-segmentation \
--segmentation-algorithm /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/E003-VPT014.json \
--input-images="/scratch/{user}/2024-01-22_MERFISH/E003-original/region_0/images/mosaic_(?P<stain>[\w|-]+)_z(?P<z>[0-9]+).tif" \
--input-micron-to-mosaic /scratch/{user}/2024-01-22_MERFISH/E003-original/region_0/images/micron_to_mosaic_pixel_transform.csv \
--output-path /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs \
--tile-size 2400 \
--tile-overlap 200
