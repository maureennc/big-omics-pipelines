#!/bin/bash

#SBATCH -n 1
#SBATCH -t 1:00:00
#SBATCH -p standard
#SBATCH -A {account}
#SBATCH --mem=4GB

# Purpose: Prepare Nanostring CosMx .tif directories/files for 'sq.read_nanostring'

# Make sure 'CellLabels' and 'CompartmentLabels' directories have been created

# Define the parent directory containing FOV subdirectories
parent_dir="/path/to/data/directory"

# Define the target directories for CellLabels and CompartmentLabels
cell_labels_dir="${parent_dir}/CellLabels"
compartment_labels_dir="${parent_dir}/CompartmentLabels"

# Loop through each FOV directory
for i in {1..130}
do
    # Format the directory name with leading zeros
    fov_dir="${parent_dir}/FOV$(printf "%03d" $i)"

    # Format the file number with leading zeros
    file_num=$(printf "%03d" $i)

    # Copy the CellLabels file
    cp "${fov_dir}/CellLabels_F${file_num}.tif" "${cell_labels_dir}/"

    # Copy the CompartmentLabels file
    cp "${fov_dir}/CompartmentLabels_F${file_num}.tif" "${compartment_labels_dir}/"
done