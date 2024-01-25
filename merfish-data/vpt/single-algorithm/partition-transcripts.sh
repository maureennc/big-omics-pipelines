#!/bin/bash

#SBATCH -A {account}
#SBATCH -p standard
#SBATCH -n 6
#SBATCH -t 02:00:00
#SBATCH --mem=40GB

module purge
module load anaconda

source activate vpt_env

vpt --verbose partition-transcripts \
--input-boundaries /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/cellpose_micron_space.parquet \
--input-transcripts /scratch/{user}/2024-01-22_MERFISH/E003-original/region_0/detected_transcripts.csv \
--output-entity-by-gene /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/E003-VPT014_cell_by_gene.csv \
--output-transcripts /scratch/{user}/2024-01-22_MERFISH/E003-VPT014/analysis_outputs/E003-VPT014_detected_transcripts.csv
