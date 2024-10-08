#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 14:14:27 2024

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
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, 'mg-sting-ko-pp-concat-hvg.h5ad'))

################################################################################################################################

# SET UP AND TRAIN MODEL A

## Model A
scvi.model.SCVI.setup_anndata(
    adata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal'],
)

model_A = scvi.model.SCVI(adata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_A.train()

################################################################################################################################

# SAVE / IMPORT MODEL

scvi_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/scvi'

## Save 
model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
#model_A.save(model_A_dir)


## Import model
scvi_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/scvi'

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
adata.layers['scvi_normalized'] = model_A.get_normalized_expression(library_size = 1e4)
adata.layers

################################################################################################################################

# INITIAL CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(adata, min_dist = 0.3)
sc.tl.leiden(adata, key_added='leiden_scVI', resolution=1)

sc.pl.umap(adata, color = ['group'])
sc.pl.umap(adata, color = ['leiden_scVI', 'Hexb'], legend_loc = 'on data')


################################################################################################################################

# GET CLUSTER GENE MARKERS

sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 30, fontsize = 6)
sc.pl.umap(adata, color = ['leiden_scVI', 'Olig1', 'Aqp4', 'Atp1a2', 'Gfap', 'Hexb', 'Rbfox3', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Cdh5', 'Rgs5', 'Tie1', 'Pdgfra'], legend_loc = 'on data')


# CLUSTER MARKERS
sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 30)
markers = sc.get.rank_genes_groups_df(adata, None)
#markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] # Keep all results
markers


################################################################################################################################

# CLUSTER VALIDATION

## Doublet score by cluster
def print_average_doublet_scores(adata, cluster_col='leiden_scVI', doublet_score_col='doublet_scores'):
    average_doublet_scores = adata.obs.groupby(cluster_col)[doublet_score_col].mean()
    for cluster, score in average_doublet_scores.items():
        print(f"Average doublet score for cluster {cluster}: {score}")

print_average_doublet_scores(adata)


genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 
         'Cx3cr1', 'Hexb', 'Mrc1',
         'Gfap', 'Aldh1l1', 'Aqp4', 'Atp1a2',
         'Pecam1', 'Tie1', 'Cldn5',
         'Pdgfrb', 'Rgs5',
         'Pdgfra', 'Cspg4', 'Sox10',
         'Mbp', 'Plp1', 'Mog'] 

sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)

################################################################################################################################


# FILTER LOW-QUALITY CLUSTERS FOR ANALYSIS

clusters_to_exclude = ['25', '28', '29', '36', '37'] # Intermediate gene expression + high average cluster doublet score
print("Clusters to exclude:", clusters_to_exclude)

### Mask and filter low-quality/non-robust clusters
mask = ~adata.obs['leiden_scVI'].isin(clusters_to_exclude)
mdata = adata[mask].copy()

print("Distribution of cells in clusters after exclusion:")
print(mdata.obs.leiden_scVI.value_counts())

################################################################################################################################

# RETRAIN MODEL ON CLEANED DATA

## Model B
scvi.model.SCVI.setup_anndata(
    mdata,
    layer='counts',
    categorical_covariate_keys=['group'],
    continuous_covariate_keys=['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal'],
)

model_B = scvi.model.SCVI(mdata)

scvi.train.Trainer(accelerator='cpu', devices=1)
model_B.train()

################################################################################################################################

# SAVING & IMPORTING

## Save model
scvi_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/scvi'
model_dir = os.path.join(scvi_dir, 'model_B') 
print(model_dir)

#model_B.save(model_dir)

## Import model
scvi_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/scvi'
model_dir = os.path.join(scvi_dir, 'model_B')
print(model_dir)

model_B = scvi.model.SCVI.load(model_dir, adata=mdata)
model_B

################################################################################################################################

# ADD LATENT KEYS

## Check for scVI entries in obsm
mdata.obsm

## add scvi latent key to obsm
SCVI_LATENT_KEY = "X_scVI"

latent = model_B.get_latent_representation()
mdata.obsm[SCVI_LATENT_KEY] = latent
latent.shape

## Add scvi normalized counts layer
mdata.layers['scvi_normalized'] = model_B.get_normalized_expression(library_size = 1e4)
mdata.layers


################################################################################################################################


# EVALUATE TRAINED MODEL B

## Extract dictionary
training_history = model_B.history
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

# SECONDARY CLUSTERING

sc.pp.neighbors(mdata, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(mdata, min_dist = 0.3)
sc.tl.leiden(mdata, key_added='leiden_scVI', resolution=1)

sc.pl.umap(mdata, color = ['leiden_scVI'], legend_loc = 'on data')


genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 
         'Cx3cr1', 'Hexb',
         'Gfap', 'Aldh1l1', 'Aqp4', 'Atp1a2',
         'Pecam1', 'Tie1', 'Cldn5',
         'Pdgfrb', 'Rgs5',
         'Pdgfra', 'Cspg4', 'Sox10',
         'Mbp', 'Plp1', 'Mog'] 

sc.tl.dendrogram(mdata, groupby = 'leiden_scVI')
sc.pl.dotplot(mdata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)

################################################################################################################################

# GET CLUSTER GENE MARKERS

sc.tl.rank_genes_groups(mdata, groupby='leiden_scVI', method='wilcoxon')
sc.pl.rank_genes_groups(mdata, n_genes = 30, fontsize = 6)
sc.pl.umap(mdata, color = ['leiden_scVI', 'Olig1', 'Aqp4', 'Atp1a2', 'Gfap', 'Hexb', 'Rbfox3', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 'Cdh5', 'Rgs5', 'Tie1', 'Pdgfra'], legend_loc = 'on data')


# CLUSTER MARKERS
sc.tl.rank_genes_groups(mdata, groupby='leiden_scVI', method='wilcoxon', layer = 'scvi_normalized', use_raw = False)
sc.pl.rank_genes_groups(mdata, n_genes = 30)
markers = sc.get.rank_genes_groups_df(mdata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] # Keep all results
markers

################################################################################################################################

# ANNOTATIONS

cell_type= { 
"0": "Excitatory neuron",
"1": "Astrocyte",
"2": "Excitatory neuron",
"3": "Oligodendrocyte",
"4": "Excitatory neuron",
"5": "Inhibitory neuron",
"6": "Inhibitory neuron",
"7": "Excitatory neuron",
"8": "Microglia",
"9": "OPC",
"10": "Inhibitory neuron",
"11": "Inhibitory neuron", 
"12": "Excitatory neuron", 
"13": "Pericyte", 
"14": "Excitatory neuron", 
"15": "Inhibitory neuron",
"16": "Oligodendrocyte", 
"17": "Endothelial cell",
"18": "Inhibitory neuron", 
"19": "Excitatory neuron",
"20": "Pericyte", 
"21": "Astrocyte", 
"22": "Excitatory neuron", 
"23": "Excitatory neuron", 
"24": "Excitatory neuron",
"25": "Inhibitory neuron",
"26": "Excitatory neuron", 
"27": "Excitatory neuron", 
"28": "Excitatory neuron", 
"29": "Astrocyte", 
"30": "Oligodendrocyte",
"31": "Excitatory neuron", 
"32": "Excitatory neuron",
"33": "Inhibitory neuron",
"34": "Inhibitory neuron",
"35": "Excitatory neuron"
}

mdata.obs['cell_type'] = mdata.obs.leiden_scVI.map(cell_type)


cluster= { 
"0": "Excitatory neuron 1",
"1": "Astrocyte 1",
"2": "Excitatory neuron 2",
"3": "Oligodendrocyte 1",
"4": "Excitatory neuron 3",
"5": "Inhibitory neuron 1",
"6": "Inhibitory neuron 2",
"7": "Excitatory neuron 4",
"8": "Microglia",
"9": "OPC",
"10": "Inhibitory neuron 3",
"11": "Inhibitory neuron 4", 
"12": "Excitatory neuron 5", 
"13": "Pericyte 1", 
"14": "Excitatory neuron 6", 
"15": "Inhibitory neuron 5",
"16": "Oligodendrocyte 2", 
"17": "Endothelial cell",
"18": "Inhibitory neuron 6", 
"19": "Excitatory neuron 7",
"20": "Pericyte 2", 
"21": "Astrocyte 2", 
"22": "Excitatory neuron 8", 
"23": "Excitatory neuron 9", 
"24": "Excitatory neuron 10",
"25": "Inhibitory neuron 7",
"26": "Excitatory neuron 11", 
"27": "Excitatory neuron 12", 
"28": "Excitatory neuron 13", 
"29": "Astrocyte 3", 
"30": "Oligodendrocyte 3",
"31": "Excitatory neuron 14", 
"32": "Excitatory neuron 15",
"33": "Inhibitory neuron 8",
"34": "Inhibitory neuron 9",
"35": "Excitatory neuron 16"
}
mdata.obs['cluster'] = mdata.obs.leiden_scVI.map(cluster)


sc.pl.umap(mdata, color = ['leiden_scVI', 'cluster'])
sc.pl.umap(mdata, color = ['cell_type'])

sc.pl.umap(mdata, color = genes, legend_loc = 'on data')


## Doublet score by cluster
def print_average_doublet_scores(mdata, cluster_col='cluster', doublet_score_col='doublet_scores'):
    average_doublet_scores = mdata.obs.groupby(cluster_col)[doublet_score_col].mean()
    for cluster, score in average_doublet_scores.items():
        print(f"Average doublet score for cluster {cluster}: {score}")

print_average_doublet_scores(mdata)
## Flag Inhibitory neuron 9


genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 
         'Cx3cr1', 'Hexb', 'Mrc1',
         'Gfap', 'Aldh1l1', 'Aqp4', 'Atp1a2',
         'Pecam1', 'Tie1', 'Cldn5',
         'Pdgfrb', 'Rgs5',
         'Pdgfra', 'Cspg4', 'Sox10',
         'Mbp', 'Plp1', 'Mog'] 


sc.pl.dotplot(mdata, genes, groupby = 'cluster', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(mdata, genes, groupby = 'cluster', standard_scale = 'var', use_raw = True, dendrogram = True)

################################################################################################################################

# CELL COMPOSITION ANALYSIS
df = pd.DataFrame(mdata.obs)

df['Cluster'] = df['cluster'].astype('category')
df['Group'] = df['dataset'].astype('category')

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
g._legend.set_title('Dataset') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()


################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/h5ad"

adata.X = csr_matrix(adata.X)
mdata.X = csr_matrix(mdata.X)

adata.write_h5ad(os.path.join(save_dir, 'mg-model_A-trained.h5ad'))
mdata.write_h5ad(os.path.join(save_dir, 'mg-model_B-cleaned-annotated.h5ad'))

#test = sc.read_h5ad(os.path.join(save_dir, 'mg-model_A-annotated.h5ad'))
#test = sc.read_h5ad(os.path.join(save_dir, 'mg-model_B-cleaned-annotated.h5ad'))

