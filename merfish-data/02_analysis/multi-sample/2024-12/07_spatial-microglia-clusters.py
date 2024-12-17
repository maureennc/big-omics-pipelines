#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 10:16:01 2024

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
import random

###############################################################################

# SETTINGS


## Matplotlib
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 500
random.seed(0)

###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))

###############################################################################

# SUBCLUSTER IMMUNE CELLS

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0)

groups = ['Microglia', 'Macrophage', 'CD4+ T cell', 'CD8+ T cell']
adata = adata[adata.obs['cell_type'].isin(groups)].copy()

sc.tl.leiden(adata, resolution=0.5, restrict_to=('cell_type', ['Microglia'])) #0.5
sc.tl.umap(adata)
sc.pl.umap(adata, color = 'leiden_R')

leiden_mapping = {
    'Macrophage': 'Macrophage',
    'CD4+ T cell': 'CD4+ T cell',
    'CD8+ T cell': 'CD8+ T cell',
    'Microglia,0': 'Mg-6',
    'Microglia,1': 'Mg-5',
    'Microglia,2': 'Mg-0',
    'Microglia,3': 'Mg-1',
    'Microglia,4': 'Mg-2',
    'Microglia,5': 'Mg-4',
    'Microglia,6': 'Mg-3'
}

adata.obs['recluster'] = adata.obs['leiden_R'].map(leiden_mapping)
sc.pl.umap(adata, color = 'recluster')
sc.pl.umap(adata, color = 'condition')



#6,5,0,1,2,4,3

###############################################################################

# SPLIT SAMPLES

E003 = adata[adata.obs['sample'] == 'E003-infected'].copy()
E007 = adata[adata.obs['sample'] == 'E007-infected'].copy()

E009 = adata[adata.obs['sample'] == 'E009-naive'].copy()

###############################################################################

# INFILTRATING SPATIAL SCATTERPLOTS

plt.rcParams['figure.figsize'] = (3,3)

## Foci
sq.pl.spatial_scatter(E003, color='inflammatory_foci', shape=None, groups = None, size = 1, img = True)

## Cells - Infected samples
sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['CD4+ T cell'], size = 1, img = True)

sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['CD8+ T cell'], size = 1, img = True)

sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['Microglia'], size = 1, img = True)

sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['Macrophage'], size = 1, img = True)

sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['CD4+ T cell', 'CD8+ T cell'], size = 1, img = True)


## Cells - naive samples
sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['CD4+ T cell'], size = 1, img = True)

sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['CD8+ T cell'], size = 1, img = True)

sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['Microglia'], size = 1, img = True)

sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['Macrophage'], size = 1, img = True)


sq.pl.spatial_scatter(E009, color='cell_type', shape=None, groups = ['CD4+ T cell', 'CD8+ T cell'], size = 1, img = True)

###############################################################################

# MICROGLIA SPATIAL SCATTERPLOTS

sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['CD4+ T cell'], size = 1, img = True)

## All microglia
sq.pl.spatial_scatter(E003, shape=None, color=["cell_type"], groups = ['Microglia'],  wspace=0.4, size = 1)


sq.pl.spatial_scatter(E003, shape=None, color=["recluster"], groups = ['Mg-1', 'CD4+ T cell'],  wspace=0.4, size =1)
sq.pl.spatial_scatter(E003, shape=None, color=["recluster"], groups = ['Mg-1'],  wspace=0.4, size = 1)



sq.pl.spatial_scatter(E003, shape=None, color=["recluster"], groups = ['Mg-2'],  wspace=0.4, size = 1)


sq.pl.spatial_scatter(E003, shape=None, color=["recluster"], groups = ['Mg-3'],  wspace=0.4, size = 1)

sq.pl.spatial_scatter(E003, shape=None, color=["recluster"], groups = ['Mg-4'],  wspace=0.4, size = 1)

sq.pl.spatial_scatter(E003, shape=None, color=["recluster"], groups = ['Mg-5'],  wspace=0.4, size = 1)

sq.pl.spatial_scatter(E003, shape=None, color=["recluster"], groups = ['Mg-6'],  wspace=0.4, size = 1)



###############################################################################

# CLUSTER MARKERS

sc.tl.rank_genes_groups(bdata, groupby='recluster')
#groups=['Mg-0', 'Mg-1', 'Mg-2', 'Mg-3', 'Mg-4', 'Mg-5', 'Mg-6']

sc.pl.rank_genes_groups(adata, groups = None)


top_genes = set()  # Using a set for uniqueness
for cluster in ['Mg-0', 'Mg-1', 'Mg-2', 'Mg-3', 'Mg-4', 'Mg-5', 'Mg-6']:
    top_genes.update(adata.uns['rank_genes_groups']['names'][cluster][:5])  # Update the set with top genes

# Convert the set back to a list
top_genes = list(top_genes)


sc.pl.heatmap(bdata, 
               var_names=top_genes, 
               groupby='recluster', 
               cmap='viridis', 
               dendrogram=True)


###############################################################################

# GENE SCORING - HOMEOSTATIC AND DAM

bdata = adata[adata.obs['cell_type'] == 'Microglia'].copy()


## Homeostatic
homeostatic_genes = ['P2ry12', 'Cx3cr1', 'Sall1', 'Sall3', 'Mertk', 'Tgfbr1', 'Tgfbr2', 'Csf1r', 'Tmem119', 'Ccr5']

sc.pl.matrixplot(bdata, var_names = homeostatic_genes, groupby = 'recluster', standard_scale = 'var', dendrogram = True)

## Dam genes
dam_genes = ['Itgax', 'Axl', 'Cybb', 'Gpnmb', 'Clec7a', 'Csf1', 'Ccl2', 'Cxcl10', 'Cxcl16', 'Arg1']
sc.pl.matrixplot(bdata, var_names = dam_genes, groupby = 'recluster', standard_scale = 'var', dendrogram = True)

## All genes
butovsky_genes = homeostatic_genes + dam_genes
sc.pl.matrixplot(bdata, var_names = butovsky_genes, groupby = 'recluster', standard_scale = 'var', dendrogram = True)


## Scoring

bdata.obs.drop(columns=['DAM score'], inplace=True, errors='ignore')
bdata.obs.drop(columns=['Homeostatic score'], inplace=True, errors='ignore')

sc.tl.dendrogram(bdata, groupby = 'leiden_R')

sc.tl.score_genes(bdata, gene_list = dam_genes)
bdata.obs.rename(columns={'score': 'DAM'}, inplace=True)

sc.tl.score_genes(bdata, gene_list = homeostatic_genes)
bdata.obs.rename(columns={'score': 'Homeostatic'}, inplace=True)


## Visualize scoring Homeostatic vs. DAM scores

### Set order
mg_order = ['Mg-0', 'Mg-1', 'Mg-2', 'Mg-3', 'Mg-4', 'Mg-5', 'Mg-6'] 
bdata.obs['recluster'] = bdata.obs['recluster'].cat.reorder_categories(mg_order)


## Plot
sc.pl.matrixplot(bdata, var_names = ['Homeostatic score', 'DAM score'], groupby = 'recluster', standard_scale = 'var', dendrogram = False, colorbar_title = 'Z-scaled scores', cmap = 'rocket', swap_axes = True)


###############################################################################

# BAR CHART OF BREAKDOWN - IMMUNE CELL COMPARTMENT COMPOSITION

## All immune cells
plt.rcParams['figure.dpi'] = 500

order = ['Mg-0', 'Mg-1', 'Mg-2', 'Mg-3', 'Mg-4', 'Mg-5', 'Mg-6', 'Macrophage', 'CD4+ T cell', 'CD8+ T cell' ]

count_data = (
    adata.obs.groupby(['recluster', 'condition'])
    .size()
    .reset_index(name='count')
)

# Plotting
hue_order = ['Naive', 'Infected']

plt.figure(figsize=(3, 3))
sns.barplot(data=count_data, x='recluster', y='count', hue='condition', palette='rocket', order = order, hue_order = hue_order)
plt.title('')
plt.xlabel('')
plt.ylabel('Count')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()

## Microglia only

order = ['Mg-0', 'Mg-1', 'Mg-2', 'Mg-3', 'Mg-4', 'Mg-5', 'Mg-6']

plt.figure(figsize=(3, 3))
sns.barplot(data=count_data, x='recluster', y='count', hue='condition', palette='rocket', order = order, hue_order = hue_order)
plt.title('')
plt.xlabel('')
plt.ylabel('Count')
plt.xticks(rotation=90)
plt.tight_layout()
plt.show()


###############################################################################

# SPATIAL PROCESSING

sq.gr.spatial_neighbors(E003)
sq.gr.nhood_enrichment(E003, cluster_key='recluster')
sq.pl.nhood_enrichment(E003, cluster_key='recluster', figsize=(2, 2), palette = 'rocket', cmap = 'rocket')


sq.gr.spatial_neighbors(E009)
sq.gr.nhood_enrichment(E009, cluster_key='recluster')
sq.pl.nhood_enrichment(E009, cluster_key='recluster', figsize=(2, 2))


###############################################################################

# EXPORT CHECKPOINT

save_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/spatial'

adata.write_h5ad(os.path.join(save_dir, 'adata-spatial.h5ad')) # immune
bdata.write_h5ad(os.path.join(save_dir, 'bdata-spatial.h5ad')) # microglia
E003.write_h5ad(os.path.join(save_dir, 'E003-spatial.h5ad'))
E009.write_h5ad(os.path.join(save_dir, 'E009-spatial.h5ad'))

###############################################################################

# IMPORT



###############################################################################

# VISUALIZATION

## Cells - Infected samples
sq.pl.spatial_scatter(E003, color='recluster', shape=None, groups = ['CD4+ T cell'], size = 1, img = True)

sq.pl.spatial_scatter(E003, color='recluster', shape=None, groups = ['CD8+ T cell'], size = 1, img = True)

sq.pl.spatial_scatter(E003, color='recluster', shape=None, groups = ['Macrophage'], size = 1, img = True)

sq.pl.spatial_scatter(E003, color='recluster', shape=None, groups = ['Mg-1', 'Mg-2', 'Mg-3', 'Mg-4', 'Mg-5', 'Mg-6'], size = 1, img = True)



###############################################################################

# QC VIOLIN PLOTS

## Set order
order = ['Naive', 'Infected'] 
adata.obs['condition'] = adata.obs['condition'].cat.reorder_categories(order)

## Plot
sc.pl.violin(adata, keys = ['total_counts'], groupby = 'sample', rotation = 90, palette = 'rocket', ylabel = 'Total counts per cell')

sc.pl.violin(adata, keys = ['n_genes_by_counts'], groupby = 'sample', rotation = 90, palette = 'rocket', ylabel = 'Detected genes per cell')

sc.pl.violin(adata, keys = ['volume'], groupby = 'sample', rotation = 90, palette = 'rocket', ylabel = 'Volume')


sc.pl.violin(adata, keys = ['total_counts'], groupby = 'condition', rotation = 90, palette = 'rocket')
sc.pl.violin(adata, keys = ['n_genes_by_counts'], groupby = 'condition', rotation = 90, palette = 'rocket')

adata.obs['n_genes_by_counts'].mean()

###############################################################################
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

sq.pl.spatial_scatter(E003, color='cell_type', shape=None, groups = ['Macrophage'], size = 1, img = True)


###############################################################################




 
#E003.obs['cropped'] = E003.obs['cropped'].astype(bool)
#E003 = E003[E003.obs['cropped'] == True].copy()