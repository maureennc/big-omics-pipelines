#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Calculate & Visualize MERFISH UMAP Embeddings from AnnData Layers

Created on Sat Dec 30 18:09:41 2023
@author: maureen
"""

##################################################################

# BEGIN

## This script calculates UMAP embeddings from different AnnData layers in a list of AnnData objects and plots them, colored by doublet score. The code is for calculating embeddings prior to running PCA.

##################################################################


import pandas as pd
import umap
import seaborn as sns
import matplotlib.pyplot as plt


##################################################################

# IMPORT DATA

## Starting point is filtered, normalized, scaled data with doublet excluded (checkpoint_A)

import os
os.chdir("/Users/maureen/Documents/data-experiments/merfish/h5ad-files/checkpoint_A/")

import scanpy as sc

adata1 = sc.read_h5ad('./adata1_post_scaling.h5ad')
adata2 = sc.read_h5ad('./adata2_post_scaling.h5ad')
adata3 = sc.read_h5ad('./adata3_post_scaling.h5ad')
adata4 = sc.read_h5ad('./adata4_post_scaling.h5ad')
adata5 = sc.read_h5ad('./adata5_post_scaling.h5ad')

adata_list = [adata1, adata2, adata3, adata4, adata5]

##################################################################

# CALCULATE EMBEDDINGS

def compute_umap_and_collect(adata, adata_index, layer_name):
    # Ensure the layer exists
    if layer_name not in adata.layers:
        raise ValueError(f"Layer {layer_name} not found in AnnData object.")

    # Use the data from the specified layer for UMAP computation
    data = adata.layers[layer_name]

    # Compute UMAP embeddings
    reducer = umap.UMAP(random_state=42, n_components=2)
    embedding = reducer.fit_transform(data)

    # Combine UMAP coordinates with cell barcodes for mapping
    umap_df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'], index=adata.obs.index)
    umap_df['Sample_Index'] = adata_index
    umap_df['Layer'] = layer_name
    umap_df['Doublet_Scores'] = adata.obs['doublet_scores'].values

    return umap_df.reset_index().rename(columns={'index': 'Cell_Barcode'})

def create_master_umap_df(adata_list, layer_names):
    all_umap_dfs = []

    for i, adata in enumerate(adata_list):
        for layer_name in layer_names:
            umap_df = compute_umap_and_collect(adata, i, layer_name)
            all_umap_dfs.append(umap_df)

    # Concatenate all DataFrames into a master DataFrame
    master_umap_df = pd.concat(all_umap_dfs, ignore_index=True)
    return master_umap_df

# Define the layers for which you want to create UMAP embeddings
layer_names = ['raw_filtered', 'normalized', 'logarithmized', 'scaled']

master_umap_df = create_master_umap_df(adata_list, layer_names)


##################################################################

# DATA VISUALIZATION


## Plot UMAP Results for each layer, colored by doublet score

def plot_umap_for_sample_and_layer(master_df, sample_index, layer_name, doublet_scores_column='Doublet_Scores'):
    # Filter the DataFrame for the specified layer and sample
    sample_layer_df = master_df[(master_df['Layer'] == layer_name) & (master_df['Sample_Index'] == sample_index)]

    # Plotting
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='UMAP1', y='UMAP2', data=sample_layer_df, hue=doublet_scores_column, palette='coolwarm', alpha=0.6)
    plt.title(f'UMAP for Sample {sample_index + 1}, {layer_name} layer')
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    plt.grid(False)
    plt.show()

def plot_all_samples_and_layers(master_df, layer_names, sample_count):
    for sample_index in range(sample_count):
        for layer_name in layer_names:
            plot_umap_for_sample_and_layer(master_df, sample_index, layer_name)

# Define the layers for which you want to create UMAP embeddings
layer_names = ['raw_filtered', 'normalized', 'logarithmized', 'scaled']

# Set number of samples
sample_count = 5

# Execute function
plot_all_samples_and_layers(master_umap_df, layer_names, sample_count)

##################################################################

# END

##################################################################