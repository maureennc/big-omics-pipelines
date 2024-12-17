#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 09:11:00 2024

@author: maureen
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
import squidpy as sq

###############################################################################

# SETTINGS


## Matplotlib
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 500

###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))


E003 = adata[adata.obs['sample'] == 'E003-infected'].copy()

#E003.obs['cropped'] = E003.obs['cropped'].astype(bool)
#E003 = E003[E003.obs['cropped'] == True].copy()

E009 = adata[adata.obs['sample'] == 'E009-naive'].copy()

###############################################################################

# SPATIAL PROCESSING

sq.gr.spatial_neighbors(E003)
sq.gr.nhood_enrichment(E003, cluster_key='cell_type')
sq.pl.nhood_enrichment(E003, cluster_key='cell_type', figsize=(3, 3))


sq.gr.spatial_neighbors(E009)
sq.gr.nhood_enrichment(E009, cluster_key='cell_type')
sq.pl.nhood_enrichment(E009, cluster_key='cell_type', figsize=(3, 3))


###############################################################################

# SPATIAL SCATTERPLOTS

import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

# Define your cell type colors
cell_type_colors = {
    'Excitatory neuron': '#FF5733',
    'Inhibitory neuron': '#C70039',
    'Oligodendrocyte': '#FFC300',
    'Astrocyte': '#DAF7A6',
    'Endothelial cell': '#900C3F',
    'Microglia': '#581845',
    'Pericyte': '#FF8D1B',
    'Macrophage': '#FF8D1B',
    'CD4+ T cell': '#5DADE2',
    'CD8+ T cell': '#3498DB',
    'OPC': '#8E44AD',
    'Choroid plexus': '#FFCE54'
}

# Create a ListedColormap
categories = list(cell_type_colors.keys())
colors = [cell_type_colors[cat] for cat in categories]
cmap = ListedColormap(colors)


## Neurons
sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['Excitatory neuron', 'Inhibitory neuron'])

sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['Excitatory neuron', 'Inhibitory neuron'])

## Glia
sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['Microglia', 'Oligodendrocyte', 'OPC', 'Astrocyte'])

sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['Microglia', 'Oligodendrocyte', 'OPC', 'Astrocyte'])


## Immune cells
sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['Macrophage', 'CD4+ T cell', 'CD8+ T cell'], size = 1, img = True)

sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['Macrophage', 'CD4+ T cell', 'CD8+ T cell'], size = 1, img = True)


###############################################################################
