#!/bin/bash

#SBATCH -A {account}
#SBATCH -p standard
#SBATCH -n 3
#SBATCH -t 03:00:00
#SBATCH --mem=60GB

module purge
module load anaconda

source activate vpt_env

vpt --verbose --processes 2 update-vzg \
--input-vzg /scratch/{user}/2024-01-22_MERFISH/E003-original/region_0/202309081147_Brain-28DaysPostInfection-090823-LS_VMSC04001_region_0.vzg \
--input-boundaries /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/cellpose_micron_space.parquet \
--input-entity-by-gene /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/E003-VPT014_cell_by_gene.csv \
--input-metadata /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/E003-VPT014_cell_metadata.csv \
--output-vzg /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/E003-VPT014.vzg
