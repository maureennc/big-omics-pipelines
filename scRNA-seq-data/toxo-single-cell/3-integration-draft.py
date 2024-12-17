#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 18:41:31 2024

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

import random
import torch
import scvi

print(sns.__version__)
print(pd.__version__)
print(np.__version__)
print(sc.__version__)
print(scvi.__version__)

################################################################################################################################

# SETTINGS

## Random seed
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 0

## Matplotlib
%matplotlib qt5
plt.rcParams['font.family'] = 'Arial'

## Scanpy
sc.set_figure_params(scanpy = True, dpi = 100, dpi_save = 400, fontsize = 14, figsize = None)

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/3/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '2-concat-pp-full.h5ad'))

################################################################################################################################

# SET UP AND TRAIN MODEL A

## Model A
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=['sample', 'group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt'], # removed pct_counts_ribosomal as covariate for simplicity
)

model_A = scvi.model.SCVI(adata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_A.train()

################################################################################################################################

# SAVE / IMPORT MODEL

scvi_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/3/scvi'

## Save 
model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
#model_A.save(model_A_dir)


## Import model
scvi_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/3/scvi'

model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
model_A = scvi.model.SCVI.load(model_A_dir, adata=adata)
model_A

################################################################################################################################

# EVALUATE TRAINED MODEL A

## Extract dictionary
training_history = model_A.history
training_history


training_history_df = pd.DataFrame(index=training_history['kl_weight'].index)

for key, df in training_history.items():
    training_history_df = training_history_df.join(df, how='outer')

## Visualize results
training_history_df.reset_index(inplace=True)

plt.figure(figsize=(5, 20))
## ELBO
plt.subplot(3, 1, 1)
plt.plot(training_history_df['epoch'], training_history_df['elbo_train'], label='ELBO')
plt.xlabel('Epochs')
plt.ylabel('ELBO')
plt.title('ELBO over Training Epochs')
plt.legend()

## Training Loss
plt.subplot(3, 1, 2)
plt.plot(training_history_df['epoch'], training_history_df['train_loss_epoch'], label='Training Loss')
plt.xlabel('Epochs')
plt.ylabel('Training Loss')
plt.title('Training Loss over Epochs')
plt.legend()

## KL Divergence (Local)
plt.subplot(3, 1, 3)
plt.plot(training_history_df['epoch'], training_history_df['kl_local_train'], label='KL Divergence (Local)')
plt.xlabel('Epochs')
plt.ylabel('KL Divergence (Local)')
plt.title('KL Divergence over Epochs')
plt.legend()

## Adjust layout
plt.tight_layout()
plt.show()


################################################################################################################################

# EXTRACT LATENT REPRESENTATION

## Check for scVI entries in obsm
adata.obsm

## add scvi latent key to obsm
SCVI_LATENT_KEY = "X_scVI"

latent = model_A.get_latent_representation()
adata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

## Add scvi normalized counts layer
adata.layers['scvi_normalized'] = model_A.get_normalized_expression()
adata.layers

################################################################################################################################

# INITIAL CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added = 'leiden_scVI', resolution = 1)

sc.pl.umap(adata, color = ['group', 'leiden_scVI'])
sc.pl.umap(adata, color = ['leiden_scVI', 'Cd163', 'Mrc1', 'Hexb', 'P2ry12', 'Sall1', 'Csf1r', 'Aif1', 'Ccr2', 'Cx3cr1', 'Apoe', 'Nos2'])

################################################################################################################################

# GET CLUSTER GENE MARKERS

sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 30)
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] # Keep all results
markers[markers['group'] == '0'].sort_values(by='logfoldchanges', ascending=False).head(50)

################################################################################################################################

# CLUSTER VALIDATION

genes = ['Cx3cr1', 'Csf1r', 'Trem2', 'Hexb', 'Sall1', 'P2ry12', 'Adora3',  'Il1a',
         'Mrc1', 'Cd163', 'Lyve1',
         'Cd3e', 'Cd3g', 'Nkg7', 'Klre1', 'Prf1', 'Cd4', 'Cd8b1', 'Gzmb',
         'Apoe', 'Cst7', 'Axl', 
         'Mmp9', 'S100a9', 'Cd19', 'Ptprc',
         'Aqp4', 'Slc2a1', 'Rbfox3'] 

sc.tl.dendrogram(adata, groupby = 'leiden_scVI')
sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)

sc.pl.umap(adata, color = 'total_counts')
sc.pl.umap(adata, color = 'leiden_scVI', legend_loc = 'on data')

################################################################################################################################

# ANNOTATIONS

cell_type= { 
"0": "Microglia", # Naive
"1": "T cell",
"2": "Macrophage",
"3": "Microglia-Macrophage",
"4": "Macrophage",
"5": "Macrophage",
"6": "Macrophage",
"7": "NK cell",
"8": "Microglia",
"9": "Macrophage",
"10": "Macrophage",
"11": "T cell",
"12": "B cell",
"13": "Macrophage",
"14": "Neutrophil",
"15": "Unknown",
"16": "T cell",
"17": "T cell"
}

adata.obs['cell_type'] = adata.obs.leiden_scVI.map(cell_type)
sc.pl.dotplot(adata, genes, groupby = 'cell_type', standard_scale = 'var', use_raw = True, dendrogram = True)


cluster= { 
"0": "Microglia 1", # Naive
"1": "T cell 1",
"2": "Macrophage 1",
"3": "Microglia-Macrophage",
"4": "Macrophage 2",
"5": "Macrophage 3",
"6": "Macrophage 4",
"7": "NK cell",
"8": "Microglia 2",
"9": "Macrophage 5",
"10": "Macrophage 6",
"11": "T cell 2",
"12": "B cell",
"13": "Macrophage 7",
"14": "Neutrophil",
"15": "Unknown",
"16": "T cell 3",
"17": "T cell 4"
    
}

adata.obs['cluster'] = adata.obs.leiden_scVI.map(cluster)

sc.pl.umap(adata, color = 'cell_type', legend_loc = 'on data')
sc.pl.umap(adata, color = 'cluster', legend_loc = 'on data')

################################################################################################################################

# VISUALIZATION

sc.pl.dotplot(adata, genes, groupby = 'cluster', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, genes, groupby = 'cluster', standard_scale = 'var', use_raw = True, dendrogram = True)

sc.pl.umap(adata, color = ['leiden_scVI'], legend_loc = 'on data', legend_fontsize = 10)
sc.pl.umap(adata, color = ['cell_type'])
sc.pl.umap(adata, color = ['cluster'])
sc.pl.umap(adata, color = ['group'])

sc.pl.umap(adata, color = ['total_counts'])
sc.pl.umap(adata, color = ['pct_counts_mt'])
sc.pl.umap(adata, color = ['pct_counts_ribosomal'])


################################################################################################################################

# CELL COMPOSITION ANALYSIS

df = pd.DataFrame(adata.obs)

df['Cell type'] = df['cell_type'].astype('category')
df['Group'] = df['group'].astype('category')
df['Cluster'] = df['cluster'].astype('category')

## Cell types
cell_type_count = df.groupby(['Group', 'Cell type']).size().reset_index(name='Count')

g = sns.catplot(
    x='Cell type',
    y='Count',
    hue='Group',  
    data=cell_type_count,
    kind='bar',
    height=5,
    aspect=2,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((1, .9))  # (x, y) The position of the legend's bounding box
g._legend.set_title('Group') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()


## Clusters
cluster_count = df.groupby(['Group', 'Cluster']).size().reset_index(name='Count')

g = sns.catplot(
    x='Cluster',
    y='Count',
    hue='Group',  
    data=cluster_count,
    kind='bar',
    height=5,
    aspect=2,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((1, .9))  # (x, y) The position of the legend's bounding box
g._legend.set_title('Group') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()

################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/3/h5ad"

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, '3-trained-annotated.h5ad'))

################################################################################################################################
