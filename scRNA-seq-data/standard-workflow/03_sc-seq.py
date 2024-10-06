#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 20 13:23:35 2024

@author: maureen

1. Train scVI model
2. Annotate leiden clusters
3. Export cluster markers DE
4. Create and export subsets for downstream analysis

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

import random
import torch
import scvi

print(sns.__version__)
print(pd.__version__)
print(np.__version__)
print(sc.__version__)
print(scvi.__version__)

###############################################################################

# SETTINGS

## Random seed
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 1

## GPU-only
#torch.cuda.manual_seed_all(0)
#torch.backends.cudnn.deterministic = True
#torch.backends.cudnn.benchmark = False

## Matplotlib
%matplotlib qt5
plt.rcParams['font.family'] = 'Arial'

## Scanpy
#sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)


###############################################################################

# IMPORT DATA

## virtual environment = 'sc-seq'

data_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/pre-processed"

adata = sc.read_h5ad(os.path.join(data_dir, "pp-concat-hvg.h5ad"))

###############################################################################

# SET UP AND TRAIN MODEL A (scVI)

## Model A
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=['sample', 'group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'n_genes_by_counts', 'doublet_scores'],
)

model_A = scvi.model.SCVI(adata)

scvi.train.Trainer(accelerator='cpu', devices=1)
#model_A.train()

###############################################################################

# SAVE / IMPORT MODEL

scvi_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/models/scvi/'

## Save 
#model_A_dir = os.path.join(scvi_dir, 'scvi_model_A')
#print(model_A_dir)
#model_A.save(model_A_dir)


## Import model
scvi_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/models/scvi'

model_A_dir = os.path.join(scvi_dir, 'scvi_model_A')
print(model_A_dir)
model_A = scvi.model.SCVI.load(model_A_dir, adata=adata)
model_A

###############################################################################

# EXTRACT LATENT REPRESENTATION

## Check for scVI entries in obsm
adata.obsm

## add scvi latent key to obsm
SCVI_LATENT_KEY = "X_scVI"

latent = model_A.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

## Add scvi normalized counts layer
adata.layers['scvi_normalized'] = model_A.get_normalized_expression(library_size = 1e4)
adata.layers

###############################################################################

# CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(adata, min_dist = 0.3)
sc.tl.leiden(adata, key_added='leiden_scVI', resolution=0.5)

sc.pl.umap(adata, color = ['leiden_scVI'], legend_loc = 'on data')
sc.pl.umap(adata, color = ['group', 'sample'])


###############################################################################

# CLUSTER MARKERS

sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 30)
sc.pl.umap(adata, color = ['leiden_scVI'])

###############################################################################

# CELL CLASS ANNOTATION 

sc.pl.umap(adata, color = ['leiden_scVI', 
                           'Cd3e', 'Trbc1', 'Klre1', 'Cd4', 'Foxp3', 'Cd8b1', 'Mki67', # T cells
                           'Klre1', # NK cells
                           'Aif1', 'Adgre1', 'Ccr2', 'Nos2' ,# Mono-mac
                           'Hexb', 'Sall1', # Microglia
                           'S100a9', # Neutrophils
                           'Cd79a', 'Mrc1'], legend_loc = 'on data')

cell_class= { 
"0": "Macrophage",
"1": "T cell",
"2": "Macrophage",
"3": "Microglia",
"4": "Microglia",
"5": "Macrophage",
"6": "NK cell",
"7": "Macrophage",
"8": "Microglia",
"9": "Macrophage",
"10": "B cell",
"11": "Neutrophil",
"12": "T cell",
"13": "Macrophage",
"14": "Macrophage",
"15": "T cell"
}

adata.obs['cell_class'] = adata.obs.leiden_scVI.map(cell_class)
sc.pl.umap(adata, color = ['leiden_scVI', 'cell_class'])

cell_type_cluster= { 
"0": "Macrophage 1",
"1": "T cell 1",
"2": "Macrophage 2",
"3": "Microglia 1",
"4": "Microglia 3",
"5": "Macrophage 3",
"6": "NK cell",
"7": "Macrophage 4",
"8": "Microglia 2",
"9": "Macrophage 5",
"10": "B cell",
"11": "Neutrophil",
"12": "T cell 2",
"13": "Macrophage 6",
"14": "Macrophage 7",
"15": "T cell 3"
}

adata.obs['cell_type_cluster'] = adata.obs.leiden_scVI.map(cell_type_cluster)
sc.pl.umap(adata, color = ['leiden_scVI', 'cell_type_cluster'], legend_loc = 'on data')

###############################################################################

# CLUSTER DIFFERENTIAL EXPRESSION SPREADSHEETS

sc.tl.rank_genes_groups(adata, groupby='cell_type_cluster', method='wilcoxon', use_raw = True)
markers = sc.get.rank_genes_groups_df(adata, group = None)
#markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] # Keep all results
markers

## Export spreadsheet
csv_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/spreadsheets/differential-expression/all"

markers.to_csv(os.path.join(csv_dir, 'scanpy-markers-DE-log1p-wilcoxon.csv'), index=False)

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/subsets-pre-hvg"

save_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/annotated-trained"

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, 'pp-annotated-model_A.h5ad'))

###############################################################################
