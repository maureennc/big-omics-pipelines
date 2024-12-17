#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 21:24:44 2024

@author: maureen
"""

# IMPORT PACKAGES

import os
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'

import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 150, dpi_save = 400)
import squidpy as sq

###############################################################################

# READ DATASETS

## Rivanna
#data_dir = "/scratch/mnc3ra/merfish-analysis/datasets"


## Local

data_dir = "/Users/maureen/Desktop/toxo-merfish/data/input-files/"

os.chdir(data_dir)


## E003, 090823-dataset, cyto2 segmentation
dir = os.path.join(data_dir, 'E003-v2')

E003 = sq.read.vizgen(path = dir, 
                       counts_file = "E003-v2_cell_by_gene.csv",
                       meta_file = "E003-v2_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")


## E007, 012624-dataset, cyto2 segmentation
dir = os.path.join(data_dir, 'E007-v2')

E007 = sq.read.vizgen(path = dir, 
                       counts_file = "E007-v2_cell_by_gene.csv",
                       meta_file = "E007-v2_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

## E008, 021624-dataset, cyto2 segmentation
dir = os.path.join(data_dir, 'E008-v1')

E008 = sq.read.vizgen(path = dir, 
                       counts_file = "E008-v1_cell_by_gene.csv",
                       meta_file = "E008-v1_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

## E009, 042524-dataset, cyto2 segmentation
dir = os.path.join(data_dir, 'E009-v1')

E009 = sq.read.vizgen(path = dir, 
                       counts_file = "E009-v1_cell_by_gene.csv",
                       meta_file = "E009-v1_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

###############################################################################

# WRITE H5AD


h5ad_dir = "/Users/maureen/Desktop/toxo-merfish/data/h5ad/raw"
os.makedirs(h5ad_dir, exist_ok=True)

# Write the AnnData objects directly using their `write` method
E003.write(f'{h5ad_dir}/E003.h5ad')
E007.write(f'{h5ad_dir}/E007.h5ad')
E008.write(f'{h5ad_dir}/E008.h5ad')
E009.write(f'{h5ad_dir}/E009.h5ad')

