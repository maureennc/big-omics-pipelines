#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 11 09:16:11 2024

@author: maureen

Exploratory Data Analysis 1

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
scvi.settings.seed = 0

## Matplotlib
%matplotlib qt5
#plt.rcParams['font.family'] = 'Arial'

## Scanpy
#sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

###############################################################################

# SET IMPORT / EXPORT DIRECTORIES

# Data and model directories
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"
scvi_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/rivanna/scvi"

###############################################################################

# IMPORT DATA

## Load base data
base_adata = sc.read_h5ad(os.path.join(data_dir, '1-tbi-seq-hvg.h5ad'))

## Define model data based on doublet score filtering
model_data = {
    'A': base_adata,
    'B': base_adata.copy(),
    'C': base_adata[base_adata.obs['doublet_scores'] < 0.3].copy(),
    'D': base_adata[base_adata.obs['doublet_scores'] < 0.3].copy(),
    'E': base_adata[base_adata.obs['doublet_scores'] < 0.1].copy(),
    'F': base_adata[base_adata.obs['doublet_scores'] < 0.2].copy(),
    'G': base_adata[base_adata.obs['doublet_scores'] < 0.2].copy(),
    'H': base_adata[base_adata.obs['doublet_scores'] < 0.1].copy()
}

## Load models corresponding to each dataset
models = {}
for label in model_data.keys():
    model_dir = os.path.join(scvi_dir, f'model_{label}')
    models[label] = scvi.model.SCVI.load(model_dir, adata=model_data[label])

## Process each model and corresponding AnnData object
for label, model in models.items():
    adata = model_data[label]

    ## Extract and store latent representation
    SCVI_LATENT_KEY = f"X_scVI_{label}"
    latent = model.get_latent_representation()
    adata.obsm[SCVI_LATENT_KEY] = latent

    ## Add normalized expression counts
    norm_expr = model.get_normalized_expression() # removed library_size=1e4
    adata.layers['scvi_normalized'] = norm_expr


###############################################################################

# MODEL ELBO PLOTS

model_viz_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/figures/eda"

for label, model in models.items():
    ## Extract training history
    training_history = model.history
    training_history_df = pd.DataFrame(index=training_history['kl_weight'].index)

    for key, df in training_history.items():
        training_history_df = training_history_df.join(df, how='outer')

    ## Visualize results
    training_history_df.reset_index(inplace=True)

    plt.figure()

    ## ELBO
    plt.subplot(3, 1, 1)
    plt.plot(training_history_df['epoch'], training_history_df['elbo_train'], label='ELBO')
    #plt.xlabel('Epochs', fontsize=10)
    #plt.ylabel('ELBO', fontsize=10)
    plt.title(f' ', fontsize=12)
    plt.legend(fontsize=8)

    ## Training Loss
    plt.subplot(3, 1, 2)
    plt.plot(training_history_df['epoch'], training_history_df['train_loss_epoch'], label='Training Loss')
    #plt.xlabel('Epochs', fontsize=10)
    #plt.ylabel('Training Loss', fontsize=10)
    plt.title(f' ', fontsize=12)
    plt.legend(fontsize=8)

    ## KL Divergence (Local)
    plt.subplot(3, 1, 3)
    plt.plot(training_history_df['epoch'], training_history_df['kl_local_train'], label='KL Divergence')
    plt.xlabel('Epochs', fontsize=10)
    #plt.ylabel('Divergence', fontsize=10)
    plt.title(f' ', fontsize=12)
    plt.legend(fontsize=8)
    plt.tight_layout(pad=10)

    ## Save
    output_file = os.path.join(model_viz_dir, f'model_{label}_history.pdf')
    plt.savefig(output_file)
    plt.close()
    print(f"Training history plot saved to {output_file}")

    
###############################################################################

# DIMENSIONALITY REDUCTION AND CLUSTERING

for label, adata in model_data.items():
    sc.pp.neighbors(adata, use_rep=f"X_scVI_{label}", random_state=0)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, key_added='leiden_scVI', resolution=1)
    print(f"model_{label} clustering", adata.obs['leiden_scVI'].value_counts())


###############################################################################

# EXPORT ANNDATA

h5ad_export_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad/eda"

for label, adata in model_data.items():
    adata.write(os.path.join(h5ad_export_dir, f'model_{label}.h5ad'))


###############################################################################

# UMAP VISUALIZATIONS - QC METRICS

os.chdir("/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/figures/eda")

umap_params = ['group', 'pct_counts_mt', 'pct_counts_ribosomal', 'total_counts', 'n_genes_by_counts', 'leiden_scVI']

for label, adata in model_data.items():
    for param in umap_params:
        sc.pl.umap(adata, color=param, save=f'_{param}_model_{label}.pdf', show=False)
        print(f"UMAP plot saved for model {label} with color parameter {param}")

for label, adata in model_data.items():
        sc.pl.umap(adata, color='doublet_scores', vmin = 0, vmax = 0.5, save=f'_doublet_scores_model_{label}.pdf', show=False)

###############################################################################

# UMAP VISUALIZATIONS - GENES OF INTEREST

#os.chdir("/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/figures/eda")

umap_genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', 
              'Cx3cr1', 'Hexb', 'P2ry12', 'Mrc1',
              'Gfap', 'Aldh1l1', 'Aqp4', 'Atp1a2',
              'Pecam1', 'Tie1', 'Cldn5',
              'Pdgfrb', 'Rgs5',
              'Pdgfra', 'Cspg4', 'Sox10',
              'Mbp', 'Plp1', 'Mog'] 

for label, adata in model_data.items():
    for gene in umap_genes:
        sc.pl.umap(adata, color=gene, save=f'_{gene}_UMAP_model_{label}.pdf', show=False)
        print(f"UMAP plot saved for model {label} with color parameter {gene}")

###############################################################################

# GENERATE CLUSTER QC SPREADSHEET

## List to hold dataframes for concatenation
all_cluster_stats = []

## Compute and collect cluster statistics
for label, adata in model_data.items():
    ## Compute various statistics by cluster
    cluster_stats = adata.obs.groupby('leiden_scVI').agg(
        num_cells=('leiden_scVI', 'size'),
        mean_doublet_scores=('doublet_scores', 'mean'),
        median_total_counts=('total_counts', 'median'),
        median_n_genes_by_counts=('n_genes_by_counts', 'median'),
        mean_pct_counts_ribosomal=('pct_counts_ribosomal', 'mean'),
        mean_pct_counts_mt=('pct_counts_mt', 'mean')
    ).reset_index()
    
    ## Add a column for the label
    cluster_stats['model'] = label
    
    ## Append to the list
    all_cluster_stats.append(cluster_stats)

## Concatenate all the dataframes
master_cluster_stats = pd.concat(all_cluster_stats, ignore_index=True)

## Save the concatenated dataframe to a CSV file
master_output_file = os.path.join(csv_export_dir, 'master_cluster_stats.csv')
master_cluster_stats.to_csv(master_output_file, index=False)
print(f"Master cluster statistics saved to {master_output_file}")


###############################################################################