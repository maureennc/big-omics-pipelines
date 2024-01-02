#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 17:14:49 2023

@author: maureen
"""

import os
from pathlib import Path
import scanpy as sc
import squidpy as sq

#################################################################

# IMPORT DATA

## original data
vizgen_dir_diameter = Path().resolve() / "merlin_watershed" / "squidpy"


adata1 = sq.read.vizgen(path = vizgen_dir_diameter, 
                       counts_file = "cell_by_gene.csv",
                       meta_file = "cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

## cyto2
vizgen_dir_diameter = Path().resolve() / "vpt_cyto2" / "squidpy"


adata2 = sq.read.vizgen(path = vizgen_dir_diameter, 
                       counts_file = "cell_by_gene.csv",
                       meta_file = "cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

##cyto2A
vizgen_dir_diameter = Path().resolve() / "vpt_cyto2A-C_diameter" / "squidpy"
adata3 = sq.read.vizgen(path = vizgen_dir_diameter, 
                       counts_file = "2A_cell_by_gene.csv",
                       meta_file = "2A_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

## cyto2B
vizgen_dir_diameter = Path().resolve() / "vpt_cyto2A-C_diameter" / "squidpy"


adata4 = sq.read.vizgen(path = vizgen_dir_diameter, 
                       counts_file = "2B_cell_by_gene.csv",
                       meta_file = "2B_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

# cyto2C
adata5 = sq.read.vizgen(path = vizgen_dir_diameter, 
                       counts_file = "2C_cell_by_gene.csv",
                       meta_file = "2C_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")


#################################################################

# BACKUP TRANSFORMATION MATRIX

os.chdir("/Users/maureen/Documents/data-experiments/merfish/planck/")

transformation_matrix = adata1.uns['spatial']['library']['scalefactors']['transformation_matrix']
transformation_matrix.to_csv('transformation_matrix_backup.csv', index=False)



## Remove unstructured data with transformation matrix
def clear_uns_in_adata_list(adata_list):
    for adata in adata_list:
        adata.uns.clear()

adata_list = [adata1, adata2, adata3, adata4, adata5]  # Replace with your list of AnnData objects
clear_uns_in_adata_list(adata_list)


#################################################################

## Save files

### individually
sc.write('../h5ad-files/908_merlin.h5ad', adata1)
sc.write('../h5ad-files/908_cyto2.h5ad', adata2)
sc.write('../h5ad-files/908_cyto2A.h5ad', adata3)
sc.write('../h5ad-files/908_cyto2B.h5ad', adata4)
sc.write('../h5ad-files/908_cyto2C.h5ad', adata5)

#############################

### loop
for i, adata in enumerate(adata_list):
    sc.write(f'adata_{i+1}.h5ad', adata)
