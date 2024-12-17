#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 10:53:50 2024

@author: maureen
"""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
import seaborn as sns

import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 150, dpi_save = 400)
from scipy.sparse import csr_matrix
import squidpy as sq
import scrublet as scr
import anndata as ad
from copy import deepcopy


################################################################################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata = sc.read_h5ad(os.path.join(data_dir, 'pp-concat.h5ad'))

################################################################################################################################

# FILTER DOUBLETS

print(f"Dataset shape before doublet filter: {adata.shape}")

threshold = 0.3
adata = adata[adata.obs['doublet_scores'] < threshold].copy()

print(f"Dataset shape after doublet filter: {filtered_adata.shape}")

################################################################################################################################

# CALCULATE SPATIAL CENTROIDS FOR EACH CELL AND PREVENT OVERLAP

## Set order
samples = ['E008-naive', 'E003-infected', 'E007-infected']

coords = {}

for sample in samples:
    # Filter data for current sample
    sample_data = adata.obs[adata.obs['sample'] == sample]
    
    # Calculate centroids
    centroid_x = (sample_data['min_x'] + sample_data['max_x']) / 2
    centroid_y = (sample_data['min_y'] + sample_data['max_y']) / 2
    
    # Assume the origin is the overall geometric center of each sample's centroids
    x_origin = centroid_x.mean()
    y_origin = centroid_y.mean()
    
    # Transform to polar coordinates
    r = np.sqrt((centroid_x - x_origin) ** 2 + (centroid_y - y_origin) ** 2)
    theta = np.arctan2(centroid_y - y_origin, centroid_x - x_origin)
    
    # Store the results back in adata.obs
    adata.obs.loc[sample_data.index, 'centroid_r'] = r
    adata.obs.loc[sample_data.index, 'centroid_theta'] = theta


################################################################################################################################

# ASSESS COVARIATE MULTILINEARITY

obs_df = adata.obs.copy()

correlation_matrix = obs_df[['total_counts', 'n_genes_by_counts', 'centroid_theta', 'centroid_r', 'volume', 'doublet_scores', 'anisotropy', 'solidity', 'perimeter_area_ratio']].corr()

print(correlation_matrix)

plt.figure(figsize=(10, 8))
sns.heatmap(correlation_matrix, annot=True, cmap='coolwarm', fmt=".2f")
plt.title('Correlation Matrix Heatmap')
plt.grid(False)
plt.show()

################################################################################################################################

# EXPORT

save_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, 'pp-polar-coordinates.h5ad'))


################################################################################################################################