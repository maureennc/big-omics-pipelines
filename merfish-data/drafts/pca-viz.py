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

# PCA SCATTERPLOTS

## Scatterplot (Default)

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


################################

## Scatterplot (Colored by Doublet Score)


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

# KERNEL DENSITY ESTIMATE PLOTS


def plot_pca_kde_plots(adata_list, n_pcs=10):
    # Use the tab10 color palette
    palette = sns.color_palette("tab10", n_colors=len(adata_list))

    # Create a large figure to hold all subplots
    fig, axes = plt.subplots(nrows=len(adata_list), ncols=n_pcs, figsize=(n_pcs * 4, len(adata_list) * 4))

    # Iterate over each AnnData object and each PC
    for i, adata in enumerate(adata_list):
        pca_df = pd.DataFrame(adata.obsm['X_pca'], columns=[f'PC{j+1}' for j in range(n_pcs)])

        for j in range(n_pcs):
            ax = axes[i][j]
            sns.kdeplot(data=pca_df, x=f'PC{j+1}', fill=True, ax=ax, color=palette[i], alpha=0.3)
            ax.set_title(f'PC{j+1}')
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_yticks([])
            ax.grid(False)  # Turn off grid lines

    # Add a grouped legend outside the plot
    #legend_labels = [f'Sample {i+1}' for i in range(len(adata_list))]
    #plt.figlegend(handles=[plt.Rectangle((0,0),1,1, color=palette[i], alpha=0.3) for i in range(len(adata_list))],
                  #labels=legend_labels, loc='upper right', bbox_to_anchor=(1.15, 0.9), title='Samples')

    # Adjust layout and aesthetics
    plt.tight_layout()
    plt.subplots_adjust(hspace=0.5, wspace=0.3)
    plt.show()

# Example usage
plot_pca_kde_plots(adata_list, n_pcs=10)


##################################################################

# HEATMAPS FOR GENE LOADINGS INTO EACH PCA

## Global Heatmap - Top Genes
## Top 50 Gene contributions to the first 10 PCs, rank ordered by row and stratified by AnnData object
### Genes with highest contributions are at the top of the plot


def plot_global_gene_loadings_heatmap(adata, n_pcs=10, sample_index=1):
    # Run PCA
    sc.tl.pca(adata, n_comps=n_pcs)
    
    # Extract PCA loadings
    loadings = adata.varm['PCs'][:, :n_pcs]
    genes = adata.var_names
    
    # Create a DataFrame with absolute loadings
    loadings_df = pd.DataFrame(np.abs(loadings), index=genes, columns=[f'PC{i+1}' for i in range(n_pcs)])

    # Sort the DataFrame by the sum of loadings across all PCs to find the top contributing genes globally
    loadings_df['Sum_Loadings'] = loadings_df.sum(axis=1)
    sorted_loadings_df = loadings_df.sort_values(by='Sum_Loadings', ascending=False).drop('Sum_Loadings', axis=1)
    
    # Plotting as a heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(sorted_loadings_df.head(50), cmap='viridis', annot=False, fmt=".2f", yticklabels=True)
    plt.title(f'Global Top Gene Contributions across First {n_pcs} PCs - Sample {sample_index}')
    plt.xlabel('Principal Component')
    plt.ylabel('Gene')
    plt.grid(False)
    plt.show()

# Example usage
for i, adata in enumerate(adata_list):
    plot_global_gene_loadings_heatmap(adata, sample_index=i+1)


#############################

## Global Heatmap, Top 50 annotated

### Genes with highest contributions are at the top of the plot


def plot_global_highest_gene_loadings_heatmap_with_annotations(adata, n_pcs=10, top_genes=50, sample_index=1):
    # Run PCA
    sc.tl.pca(adata, n_comps=n_pcs)
    
    # Extract PCA loadings
    loadings = adata.varm['PCs'][:, :n_pcs]
    genes = adata.var_names
    
    # Create a DataFrame with loadings
    loadings_df = pd.DataFrame(loadings, index=genes, columns=[f'PC{i+1}' for i in range(n_pcs)])

    # Calculate the sum of absolute loadings across all PCs
    loadings_df['Sum_Loadings'] = loadings_df.abs().sum(axis=1)

    # Sort the DataFrame by the sum of absolute loadings in descending order
    sorted_loadings_df = loadings_df.sort_values(by='Sum_Loadings', ascending=False).drop('Sum_Loadings', axis=1)
    
    # Define custom annotation parameters
    annot_params = {"fontsize": 8}  # Adjust the fontsize as needed
    
    # Plotting as a heatmap with annotations
    plt.figure(figsize=(15, 9))
    sns.heatmap(sorted_loadings_df.head(top_genes), cmap='inferno', annot=True, fmt=".2f", yticklabels=True, annot_kws=annot_params)
    plt.title(f'Global Highest Gene Contributions across First {n_pcs} PCs - Sample {sample_index}')
    plt.xlabel('Principal Component')
    plt.ylabel('Gene')
    plt.yticks(fontsize = 8)
    plt.grid(False)
    plt.show()

# Example usage
for i, adata in enumerate(adata_list):
    plot_global_highest_gene_loadings_heatmap_with_annotations(adata, sample_index=i+1, top_genes=50)

#############################    

## Global Heatmap - Bottom Genes

### Genes with lowest contributions are at the top of the plot

def plot_global_lowest_gene_loadings_heatmap(adata, n_pcs=10, sample_index=1):
    # Run PCA
    sc.tl.pca(adata, n_comps=n_pcs)
    
    # Extract PCA loadings
    loadings = adata.varm['PCs'][:, :n_pcs]
    genes = adata.var_names
    
    # Create a DataFrame with absolute loadings
    loadings_df = pd.DataFrame(np.abs(loadings), index=genes, columns=[f'PC{i+1}' for i in range(n_pcs)])

    # Sort the DataFrame by the sum of loadings across all PCs to find the lowest contributing genes globally
    loadings_df['Sum_Loadings'] = loadings_df.sum(axis=1)
    sorted_loadings_df = loadings_df.sort_values(by='Sum_Loadings', ascending=True).drop('Sum_Loadings', axis=1)
    
    # Plotting as a heatmap
    plt.figure(figsize=(12, 8))
    sns.heatmap(sorted_loadings_df.head(50), cmap='viridis', annot=False, fmt=".2f", yticklabels=True)
    plt.title(f'Global Lowest Gene Contributions across First {n_pcs} PCs - Sample {sample_index}')
    plt.xlabel('Principal Component')
    plt.ylabel('Gene')
    plt.grid(False)
    plt.show()

# Example usage
for i, adata in enumerate(adata_list):
    plot_global_lowest_gene_loadings_heatmap(adata, sample_index=i+1)


##############################

## Global Heatmap - Bottom Genes, Annotated

### Genes with lowest contributions are at the top of the plot

def plot_global_lowest_gene_loadings_heatmap_with_annotations(adata, n_pcs=10, bottom_genes=50, sample_index=1):
    # Run PCA
    sc.tl.pca(adata, n_comps=n_pcs)
    
    # Extract PCA loadings
    loadings = adata.varm['PCs'][:, :n_pcs]
    genes = adata.var_names
    
    # Create a DataFrame with absolute loadings
    loadings_df = pd.DataFrame(loadings, index=genes, columns=[f'PC{i+1}' for i in range(n_pcs)])

    # Calculate the sum of absolute loadings across all PCs
    loadings_df['Sum_Loadings'] = loadings_df.abs().sum(axis=1)

    # Sort the DataFrame by the sum of absolute loadings in ascending order
    sorted_loadings_df = loadings_df.sort_values(by='Sum_Loadings', ascending=True).drop('Sum_Loadings', axis=1)
    
    # Define custom annotation parameters
    annot_params = {"fontsize": 8}  # Adjust the fontsize as needed
    
    # Plotting as a heatmap with annotations
    plt.figure(figsize=(15, 9))
    sns.heatmap(sorted_loadings_df.head(bottom_genes), cmap='viridis', annot=True, fmt=".2f", yticklabels=True, annot_kws=annot_params)
    plt.title(f'Global Lowest Gene Contributions across First {n_pcs} PCs - Sample {sample_index}')
    plt.xlabel('Principal Component')
    plt.ylabel('Gene')
    plt.yticks(fontsize = 8)
    plt.grid(False)
    plt.show()

# Example usage
for i, adata in enumerate(adata_list):
    plot_global_lowest_gene_loadings_heatmap_with_annotations(adata, sample_index=i+1, bottom_genes=50)


##################################################################

# Bar Graph (visualizing individual PCs)

## Top 50 genes that go into each PC

def plot_top_50_gene_loadings_bar(adata, n_pcs=10, top_genes=50, sample_index=1):
    # Run PCA
    sc.tl.pca(adata, n_comps=n_pcs)

    # Extract PCA loadings
    loadings = adata.varm['PCs'][:, :n_pcs]
    genes = adata.var_names

    # Create a DataFrame with loadings
    loadings_df = pd.DataFrame(loadings, index=genes, columns=[f'PC{i+1}' for i in range(n_pcs)])

    # Plotting a bar graph for each PC with top 50 genes
    for pc in range(n_pcs):
        # Select top 50 genes for the current PC by absolute loading
        top_genes_indices = loadings_df[f'PC{pc+1}'].abs().nlargest(top_genes).index
        top_genes_df = loadings_df.loc[top_genes_indices, [f'PC{pc+1}']]
        top_genes_df = top_genes_df.reset_index()
        top_genes_df.columns = ['Gene', 'Loading']

        # Seaborn barplot
        plt.figure(figsize=(10, 6))
        sns.barplot(y='Gene', x='Loading', data=top_genes_df.sort_values(by='Loading', ascending=True), palette='viridis')
        plt.title(f'Top 50 Gene Contributions to PC{pc+1} - Sample {sample_index}')
        plt.xlabel('Loading')
        plt.ylabel('Gene')
        plt.yticks(fontsize = 8)
        plt.grid(False)
        plt.show()

# Example usage
for i, adata in enumerate(adata_list):
    plot_top_50_gene_loadings_bar(adata, sample_index=i+1)


####################################

# Bar Graph (visualizing individual PCs)

## Bottom 50 genes that go into each PC

def plot_bottom_50_gene_loadings_bar(adata, n_pcs=10, bottom_genes=50, sample_index=1):
    # Run PCA
    sc.tl.pca(adata, n_comps=n_pcs)

    # Extract PCA loadings
    loadings = adata.varm['PCs'][:, :n_pcs]
    genes = adata.var_names

    # Create a DataFrame with loadings
    loadings_df = pd.DataFrame(loadings, index=genes, columns=[f'PC{i+1}' for i in range(n_pcs)])

    # Plotting a bar graph for each PC with bottom 50 genes
    for pc in range(n_pcs):
        # Select bottom 50 genes for the current PC by absolute loading
        bottom_genes_indices = loadings_df[f'PC{pc+1}'].abs().nsmallest(bottom_genes).index
        bottom_genes_df = loadings_df.loc[bottom_genes_indices, [f'PC{pc+1}']]
        bottom_genes_df = bottom_genes_df.reset_index()
        bottom_genes_df.columns = ['Gene', 'Loading']

        # Seaborn barplot
        plt.figure(figsize=(10, 6))
        sns.barplot(y='Gene', x='Loading', data=bottom_genes_df.sort_values(by='Loading'), palette='magma')
        plt.title(f'Bottom 50 Gene Contributions to PC{pc+1} - Sample {sample_index}')
        plt.xlabel('Loading')
        plt.ylabel('Gene')
        plt.yticks(fontsize=8)
        plt.grid(False)
        plt.show()

# Example usage
for i, adata in enumerate(adata_list):
    plot_bottom_50_gene_loadings_bar(adata, sample_index=i+1)



###########################################################################

# Bar Plot, Loadings for Top 20% of Genes

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

# PCA, Concatenated Results; Visualizes individual cells

import anndata as ad

all_data = ad.concat(adata_list, join='outer', label='batch', keys=['adata1', 'adata2', 'adata3', 'adata4', 'adata5'])

# Perform PCA on the combined dataset
sc.tl.pca(all_data)

# Visualize the results
sc.pl.pca_scatter(all_data, color='batch')



##################################################################

