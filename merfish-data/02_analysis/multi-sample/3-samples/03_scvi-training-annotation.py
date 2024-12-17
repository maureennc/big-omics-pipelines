#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 20:02:35 2024

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

################################################################################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata = sc.read_h5ad(os.path.join(data_dir, 'pp-polar-coordinates.h5ad'))

################################################################################################################################

# SET UP AND TRAIN MODEL A

## Model A
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=['sample', 'condition'],
    continuous_covariate_keys=['total_counts', 'doublet_scores', 'centroid_theta', 'centroid_r', 'anisotropy'],
)

model_A = scvi.model.SCVI(adata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_A.train()


################################################################################################################################

# SAVE / IMPORT MODEL

scvi_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/scvi'

## Save 
model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
#model_A.save(model_A_dir)


## Import model
scvi_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/scvi'

model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
model_A = scvi.model.SCVI.load(model_A_dir, adata=adata)
model_A

################################################################################################################################

# EVALUATE TRAINED MODEL B

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
adata.layers['scvi_normalized'] = model_A.get_normalized_expression(library_size = 1e4)
adata.layers

################################################################################################################################

# INITIAL CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(adata, min_dist = 0.3)
sc.tl.leiden(adata, key_added='leiden_scVI', resolution=1)

sc.pl.umap(adata, color = ['leiden_scVI'], legend_loc = 'on data')
sc.pl.umap(adata, color = ['leiden_scVI', 'condition', 'sample'])

sc.pl.umap(adata, color = ['P2ry12', 'Tnfaip2', 'Ccr2', 'Cd4', 'Cd8a', 'Rbfox3', 'Olig1', 'Pdgfra', 'Rgs5', 'Pecam1', 'Tie1', 'Aqp4', 'Slc17a6', 'Slc17a7', 'Gad2', 'Gas6'])


################################################################################################################################

# GET CLUSTER GENE MARKERS

sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon')
sc.pl.rank_genes_groups(adata, n_genes = 30, fontsize = 12)


markers = sc.get.rank_genes_groups_df(adata, None)
#markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] # Keep all results
markers


genes = ['Rbfox3', 'Slc17a6', 'Slc17a7', 'Gad2','Ache',
         'Cx3cr1', 'P2ry12', 'Ccr2', 'Tnfaip2',
         'Aldh1l1', 'Aqp4', 'Atp1a2',
         'Pecam1', 'Tie1', 'Cldn5',
         'Rgs5', 'Gas6', 'Cd4', 'Cd8a',
         'Pdgfra', 'Cspg4', 'Sox10'] 

sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)


################################################################################################################################

# CLUSTER VALIDATION

## Doublet score by cluster
def print_average_doublet_scores(adata, cluster_col='leiden_scVI', doublet_score_col='doublet_scores'):
    average_doublet_scores = adata.obs.groupby(cluster_col)[doublet_score_col].mean()
    for cluster, score in average_doublet_scores.items():
        print(f"Average doublet score for cluster {cluster}: {score}")

print_average_doublet_scores(adata)

sc.tl.dendrogram(adata, groupby = 'cluster')
sc.pl.umap(adata, color = 'leiden_scVI', legend_loc = 'on data')
sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)


################################################################################################################################

# CELL CLASS ANNOTATION 


cell_class= { 
"0": "Neuron", # Excitatory
"1": "Oligodendrocyte",
"2": "Astrocyte", 
"3": "Vascular", # Endothelial cell
"4": "Neuron", # Inhibitory 
"5": "Microglia", 
"6": "Neuron", # Excitatory
"7": "CD4+ T cell",
"8": "Macrophage",
"9": "Neuron", # Excitatory 
"10": "CD8+ T cell", 
"11": "Neuron", # Inhibitory 
"12": "Neuron", # Excitatory
"13": "Vascular", # Endothelial-Pericyte intermediate 
"14": "Neuron", # Excitatory 
"15": "Oligodendrocyte",
"16": "OPC",
"17": "Neuron", # Excitatory
"18": "Choroid plexus",
"19": "Choroid plexus",
"20": "Neuron" # Excitatory
}

adata.obs['cell_class'] = adata.obs.leiden_scVI.map(cell_class)
sc.pl.umap(adata, color = ['leiden_scVI', 'cell_class'], legend_loc = 'on data')

cell_type= { 
"0": "Excitatory neuron",
"1": "Oligodendrocyte",
"2": "Astrocyte", 
"3": "Vascular", # Endothelial cell
"4": "Inhibitory neuron", 
"5": "Microglia", 
"6": "Excitatory neuron",
"7": "CD4+ T cell",
"8": "Macrophage",
"9": "Excitatory neuron", 
"10": "CD8+ T cell", 
"11": "Inhibitory neuron", 
"12": "Excitatory neuron",
"13": "Vascular", # Endothelial-Pericyte intermediate 
"14": "Excitatory neuron", 
"15": "Oligodendrocyte",
"16": "OPC",
"17": "Excitatory neuron",
"18": "Choroid plexus",
"19": "Choroid plexus",
"20": "Excitatory neuron" # Excitatory
}

adata.obs['cell_type'] = adata.obs.leiden_scVI.map(cell_type)
sc.pl.umap(adata, color = ['cell_type'], legend_loc = 'on data')


cluster= { 
"0": "Excitatory neuron 1",
"1": "Oligodendrocyte 1",
"2": "Astrocyte", 
"3": "Vascular 1", # Endothelial cell-dominant
"4": "Inhibitory neuron 1", 
"5": "Microglia", 
"6": "Excitatory neuron 2",
"7": "CD4+ T cell",
"8": "Macrophage",
"9": "Excitatory neuron 3", 
"10": "CD8+ T cell", 
"11": "Inhibitory neuron 2", 
"12": "Excitatory neuron 4",
"13": "Vascular 2", # Pericyte-dominant 
"14": "Excitatory neuron 5", 
"15": "Oligodendrocyte 2",
"16": "OPC",
"17": "Excitatory neuron 6",
"18": "Choroid plexus",
"19": "Choroid plexus",
"20": "Excitatory neuron 7"
}

adata.obs['cluster'] = adata.obs.leiden_scVI.map(cluster)
sc.pl.umap(adata, color = ['cluster'], legend_loc = 'on data')


################################################################################################################################

# CELL COMPOSITION ANALYSIS
df = pd.DataFrame(adata.obs)

df['Cluster'] = df['cluster'].astype('category')
df['Group'] = df['sample'].astype('category')

count_df = df.groupby(['Group', 'Cluster']).size().reset_index(name='Count')

plt.close()
g = sns.catplot(
    x='Cluster',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=2,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((1.05, 1))  # (x, y) The position of the legend's bounding box
g._legend.set_title('sample') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()



df['Cell type'] = df['cell_type'].astype('category')
df['Group'] = df['sample'].astype('category')

count_df = df.groupby(['Group', 'Cell type']).size().reset_index(name='Count')

plt.close()
g = sns.catplot(
    x='Cell type',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=2,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((1.05, 1))  # (x, y) The position of the legend's bounding box
g._legend.set_title('sample') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)

plt.show()

################################################################################################################################

# EXPORT

save_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, 'scvi-trained_A-annotated.h5ad'))

################################################################################################################################
