#!/bin/bash

#SBATCH -A {account}
#SBATCH -p standard
#SBATCH -n 4
#SBATCH -t 01:00:00
#SBATCH --mem=10GB

module purge
module load anaconda

source activate vpt_env

vpt --verbose derive-entity-metadata \
--input-boundaries /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/cellpose_micron_space.parquet \
--input-entity-by-gene /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/E003-VPT014_cell_by_gene.csv \
--output-metadata /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/E003-VPT014_cell_metadata.csv
