#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 10:31:08 2024

@author: maureen
"""

import os

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 100, dpi_save = 400)

import squidpy as sq
import scrublet as scr
import umap

################################################################################################################################

# IMPORT DATA

base_dir = "/Users/maureen/Documents/projects/harris-lab/MERSCOPE/data/vpt-datasets/E003_20230908-infected/"

# List of dataset suffixes
start = 29
end = 37
dataset_suffixes = [str(i).zfill(3) for i in range(start, end + 1)]

# Initialize list
adata_list = []

for suffix in dataset_suffixes:
    # Construct the dataset directory and file names
    dataset_dir = f"E003-VPT{suffix}"
    counts_file = f"E003-VPT{suffix}_cell_by_gene.csv"
    meta_file = f"E003-VPT{suffix}_cell_metadata.csv"
    transformation_file = "micron_to_mosaic_pixel_transform.csv"

    # Construct the full paths
    counts_path = f"{base_dir}{dataset_dir}/{counts_file}"
    meta_path = f"{base_dir}{dataset_dir}/{meta_file}"
    transformation_path = f"{base_dir}{dataset_dir}/{transformation_file}"

    # Directly assign the AnnData object to a dynamically named variable
    globals()[f'adata_{suffix}'] = sq.read.vizgen(path=base_dir + dataset_dir, 
                                                  counts_file=counts_file,
                                                  meta_file=meta_file,
                                                  transformation_file=transformation_file)
    
    # Append the AnnData object to the list
    adata_list.append(globals()[f'adata_{suffix}'])


################################################################################################################################

# QC 1

def plot_qc1(adata_list):

    qc_df = pd.DataFrame()  # Initialize the master DataFrame

    for i, adata in enumerate(adata_list):
        # QC Calculations
        sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)
        
        # Selecting relevant QC metrics
        qc_metrics = ['total_counts', 'n_genes_by_counts', 'volume']
        df_qc = adata.obs[qc_metrics].copy()
        df_qc = df_qc.melt(var_name='QC Metric', value_name='Value')
        df_qc['Sample'] = f'Sample {i+1}'

        # Append to the master DataFrame
        qc_df = pd.concat([qc_df, df_qc], ignore_index=True)

        # Plotting for the current sample
        g = sns.FacetGrid(df_qc, col='QC Metric', sharex=False, sharey=False)
        g.map_dataframe(sns.histplot, x='Value')

        # Adjusting subplot titles
        for ax in g.axes.flatten():
            ax.set_title(ax.get_title(), fontsize=10)  # Set smaller font size here

        g.fig.suptitle(f'QC1 for E003-VPT{i+ 29}', fontsize=16)
        g.fig.subplots_adjust(top=0.85)
        plt.show()

    return qc_df


qc_df = plot_qc1(adata_list)
del qc_df #clear qc_df from memory


################################################################################################################################

# FILTER POOR QUALITY CELLS


def pp_and_filter(adata_list, min_volume, max_volume):
    for adata in adata_list:
        #adata.layers['raw_import'] = adata.X.copy() # Copy raw data to a separate layer
        
        # Applying initial filters
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.filter_cells(adata, min_genes=10)
        sc.pp.filter_cells(adata, min_counts=50)
        sc.pp.filter_cells(adata, max_counts=3000)
        
        # Apply volume filter using boolean indexing
        volume_mask = adata.obs['volume'].between(min_volume, max_volume)
        adata._inplace_subset_obs(volume_mask)

        # Calculate QC metrics
        sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)


min_volume = 50
max_volume = 4000
pp_and_filter(adata_list, min_volume, max_volume)



################################################################################################################################

# QC2

def plot_qc2(adata_list):
    qc_df = pd.DataFrame()  # Initialize the master DataFrame

    for i, adata in enumerate(adata_list):
        # Recalculate QC Metrics after filtering
        sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)
        
        # Select relevant QC metrics
        qc_metrics = ['total_counts', 'n_genes_by_counts', 'volume']
        df_qc = adata.obs[qc_metrics].copy()
        df_qc = df_qc.melt(var_name='QC Metric', value_name='Value')
        df_qc['Sample'] = f'Sample {i+1}'

        # Append to master DataFrame
        qc_df = pd.concat([qc_df, df_qc], ignore_index=True)

        # Plotting for the current sample
        g = sns.FacetGrid(df_qc, col='QC Metric', sharex=False, sharey=False)
        g.map_dataframe(sns.histplot, x='Value')

        # Adjusting subplot titles
        for ax in g.axes.flatten():
            ax.set_title(ax.get_title(), fontsize=10)  # Set smaller font size here

        g.fig.suptitle(f'QC2 for VPT {i+29}', fontsize=16)
        g.fig.subplots_adjust(top=0.85)
        plt.show()

    return qc_df

qc_df = plot_qc2(adata_list)
del qc_df #clear qc_df from memory


################################################################################################################################

# DOUBLET DETECTION ON FILTERED DATA


def run_scrublet(adata_list):
    # Initialize an empty list to store DataFrame rows
    scrublet_rows = []

    for i, adata in enumerate(adata_list):
        # Run Scrublet
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.20)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=10) #decreased from 30
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        # Plot Histograms
        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for VPT{i+29}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(False)
        plt.show()

        # Extract cell barcodes
        cell_barcodes = adata.obs.index

        # Store Scrublet data with cell barcodes for each AnnData object in the list
        for barcode, obs_score, sim_score, pred_doublet in zip(cell_barcodes, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets):
            scrublet_rows.append({'Sample_Index': i+1, 
                                  'Cell_Barcode': barcode,
                                  'Observed_Score': obs_score, 
                                  'Simulated_Score': sim_score, 
                                  'Predicted_Doublet': pred_doublet})

    # Create a DataFrame from the list of rows
    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df

# Run Scrublet analysis
scrublet_result = run_scrublet(adata_list)



################################################################################################################################

# View doublet scores
# FILTER DOUBLETS

## Doublet counts

# Assuming 'adata_list' is your list of AnnData objects and they all have the 'predicted_doublets' column in their .obs
start = 29
end = 37
sample_names = [f'VPT{str(i).zfill(3)}' for i in range(start, end + 1)]

# Initialize lists to store counts
true_doublets = []
false_doublets = []

# Extracting the counts of predicted doublets
for adata in adata_list:
    value_counts = adata.obs.predicted_doublets.value_counts()
    true_doublets.append(value_counts.get(True, 0))  # Get count of True, default to 0 if not found
    false_doublets.append(value_counts.get(False, 0))  # Get count of False, default to 0 if not found

# Plotting
x = range(len(adata_list))  # X-axis points
width = 0.35  # Width of the bars

fig, ax = plt.subplots()
rects1 = ax.bar(x, false_doublets, width, label='Not Doublets')
rects2 = ax.bar(x, true_doublets, width, bottom=false_doublets, label='Predicted Doublets')

# Add some text for labels, title, and custom x-axis tick labels, etc.
ax.set_ylabel('Counts')
ax.set_title('Counts by sample and doublet status')
ax.set_xticks(x)
ax.set_xticklabels(sample_names, rotation=45, ha='right')
ax.legend()
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.grid(False)
plt.show()



## Perform filtering
def filter_doublets(adata_list):
    for adata in adata_list:
        # Filter out cells where 'predicted_doublets' is True
        adata._inplace_subset_obs(adata.obs['predicted_doublets'] == False)

# Apply the doublet filtering step
filter_doublets(adata_list)


################################################################################################################################

# SAVE LAYER OF RAW, FILTERED DATA FOR EACH ANNDATA OBJECT


def save_raw_to_layer(adata_list, layer_name='raw_filtered'):
    for adata in adata_list:
        adata.layers[layer_name] = adata.X.copy()

save_raw_to_layer(adata_list)



################################################################################################################################

def preprocess_adata(adata, layers_dict):

    # Normalization
    sc.pp.normalize_total(adata, target_sum=1e4)
    adata.layers[layers_dict['normalized']] = adata.X.copy()

    # Logarithmization
    sc.pp.log1p(adata)
    adata.layers[layers_dict['logarithmized']] = adata.X.copy()

    # Scaling
    sc.pp.scale(adata, max_value=10)
    adata.layers[layers_dict['scaled']] = adata.X.copy()

def process_adata_list(adata_list):
    layers_dict = {
        'normalized': 'normalized',
        'logarithmized': 'logarithmized',
        'scaled': 'scaled'
    }

    for adata in adata_list:
        preprocess_adata(adata, layers_dict)

process_adata_list(adata_list)

################################################################################################################################

# CLUSTERING

## 10 PCs, 15 neighbors
for i, adata in enumerate(adata_list):
    sc.tl.pca(adata, n_comps = 10)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=10)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5)
    sc.pl.umap(adata, color=['leiden'])

