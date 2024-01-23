#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 13:10:00 2024

@author: maureen
"""

# SPATIAL SEGMENT

## Plot single segment

sq.pl.spatial_segment(
    adata,
    color="leiden_scVI",
    library_key="fov",
    library_id=['1'], # Toggle through FOVs
    seg_cell_id="cell_ID",
)

## Export each FOV with segmentation masks as .png
output_directory = '/Users/maureen/Documents/projects/lukens-lab/nick/nanostring-cosmx/spatial_segment_FOVs'

for fov_number in range(1, 131):
    # Call sq.pl.spatial_segment to generate the plot for the current FOV
    plot = sq.pl.spatial_segment(
        adata,
        color="leiden_scVI",
        library_key="fov",
        library_id=[str(fov_number)], 
        seg_cell_id="cell_ID",
    )

    file_name = f"fov_{fov_number}_spatial_segment.png"
    save_path = os.path.join(output_directory, file_name)

    plt.savefig(save_path)
    plt.close(plot)
