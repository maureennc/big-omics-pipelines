#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 30 22:17:41 2023

@author: maureen
"""


import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

##################################################################

# IMPORT DATA

## Starting point is data that has been filtered, normalized, logarithmized, and scaled (checkpoint_A)

import os
os.chdir("/Users/maureen/Documents/data-experiments/merfish/h5ad-files/checkpoint_A/")

adata1 = sc.read_h5ad('./adata1_post_scaling.h5ad')
adata2 = sc.read_h5ad('./adata2_post_scaling.h5ad')
adata3 = sc.read_h5ad('./adata3_post_scaling.h5ad')
adata4 = sc.read_h5ad('./adata4_post_scaling.h5ad')
adata5 = sc.read_h5ad('./adata5_post_scaling.h5ad')

adata_list = [adata1, adata2, adata3, adata4, adata5]


##################################################################

# RUN PCA

def run_pca(adata_list, n_pcs=50):
    for adata in adata_list:
        sc.pp.pca(adata, n_comps=n_pcs)
    
    return adata_list #adata.varm should be added

run_pca(adata_list) 

##################################################################

# SCREE PLOT

def plot_scree(adata, title=''):
    explained_variance_ratio = adata.uns['pca']['variance_ratio']

    plt.figure(figsize=(10, 6))
    plt.plot(range(1, len(explained_variance_ratio) + 1), explained_variance_ratio, 'o-')
    plt.title(title)
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Explained')
    plt.show()

for i, adata in enumerate(adata_list):
    plot_scree(adata, title=f'Scree Plot for Sample {i+1}')




##################################################################

# PCA Scatterplot (Default)

def plot_pca_scatter(adata, title=''):
    pca_df = pd.DataFrame(adata.obsm['X_pca'], columns=[f'PC{i+1}' for i in range(adata.obsm['X_pca'].shape[1])])
    pca_df['cell_id'] = adata.obs_names

    plt.figure(figsize=(8, 6))
    sns.scatterplot(x='PC1', y='PC2', data=pca_df, s=50, alpha=0.7)
    plt.title(title)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.grid(False)
    plt.show()

for i, adata in enumerate(adata_list):
    plot_pca_scatter(adata, title=f'PCA Scatter Plot for Sample {i+1}')


##################################################################

# PCA Scatterplot (DOUBLET SCORE)


def plot_pca_scatter(adata, title=''):
    # Create a DataFrame from the PCA results
    pca_df = pd.DataFrame(adata.obsm['X_pca'], columns=[f'PC{i+1}' for i in range(adata.obsm['X_pca'].shape[1])])
    pca_df['cell_id'] = adata.obs_names

    # Ensure that the cell_id in pca_df matches the index of adata.obs
    # Add doublet scores to the DataFrame by mapping based on cell_id
    pca_df = pca_df.set_index('cell_id')
    pca_df['doublet_scores'] = adata.obs['doublet_scores']
    
    plt.figure(figsize=(8, 6))
    # Use doublet scores for coloring the scatter plot
    scatter = plt.scatter(pca_df['PC1'], pca_df['PC2'], c=pca_df['doublet_scores'], s=50, alpha=0.7, cmap='viridis')
    plt.title(title)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.colorbar(scatter, label='Doublet Scores')  # Add a colorbar for doublet scores
    plt.grid(False)
    plt.show()

for i, adata in enumerate(adata_list):
    plot_pca_scatter(adata, title=f'PCA Scatter Plot for Sample {i+1}')


##################################################################

# PCA LOADINGS BAR PLOT
def plot_top_pca_genes_facet(adata, n_pcs=10, top_percent=0.20, sample_index=1):
    # Run PCA
    sc.tl.pca(adata, n_comps=n_pcs)

    # Number of genes to show (top 20%)
    num_genes_to_show = int(np.ceil(top_percent * adata.var.shape[0]))

    # Create a figure with subplots
    fig, axes = plt.subplots(nrows=2, ncols=5, figsize=(20, 10))  # Adjust the size as needed
    axes = axes.flatten()  # Flatten the array for easy iteration

    for pc in range(n_pcs):
        # Extract PCA loadings
        loadings = adata.varm['PCs'][:, pc]
        top_genes_indices = np.argsort(np.abs(loadings))[-num_genes_to_show:]
        top_genes = adata.var_names[top_genes_indices]
        top_loadings = loadings[top_genes_indices]

        # Create DataFrame for plotting
        loading_df = pd.DataFrame({
            'Gene': top_genes,
            'Loading': top_loadings
        }).sort_values(by='Loading', ascending=False)

        # Plotting on subplot
        ax = axes[pc]
        sns.barplot(x='Loading', y='Gene', data=loading_df, ax=ax)
        ax.set_title(f'PC{pc+1}')
        ax.set_xlabel('PCA Loading')
        ax.set_ylabel('')
        ax.tick_params(axis='y', labelsize=6)  # Adjust font size

    plt.tight_layout()
    fig.suptitle(f'Sample {sample_index} - Top 20% Gene Contributions to PCs', fontsize=16)
    plt.subplots_adjust(top=0.9)  # Adjust the top padding
    plt.show()

# Example usage
for i, adata in enumerate(adata_list):
    plot_top_pca_genes_facet(adata, sample_index=i+1)



##################################################################



# NEIGHBORS

for adata in adata_list:
    sc.pp.neighbors(adata, n_neighbors=30, use_rep='X_pca', metric='euclidean')
    sc.tl.umap(adata)  # for visualization
    sc.tl.leiden(adata)  # for clustering
    # Visualize and assess the results



##################################################################

