#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  2 12:15:55 2024

@author: maureen
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 150, dpi_save = 400)
import scrublet as scr

# environment = 'sc'

##################################################################

# IMPORT DATA


os.chdir("/Users/maureen/Documents/data-experiments/merfish/h5ad-files")

## these files were processed as in 'save-anndata-as-h5ad'
adata1 = sc.read_h5ad('908_merlin.h5ad')
adata2 = sc.read_h5ad('908_cyto2.h5ad')
adata3 = sc.read_h5ad('908_cyto2A.h5ad')
adata4 = sc.read_h5ad('908_cyto2B.h5ad')
adata5 = sc.read_h5ad('908_cyto2C.h5ad')

##################################################################

# QC 1

##################################################################

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


adata_list = [adata1, adata2, adata3, adata4, adata5]
min_volume = 50
max_volume = 4000
pp_and_filter(adata_list, min_volume, max_volume)


##################################################################

# QC2

# QC3

##################################################################

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
                                                                  n_prin_comps=30)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        # Plot Histograms
        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for AnnData Object {i+1}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
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
adata_list = [adata1, adata2, adata3, adata4, adata5]  # Replace with your list of AnnData objects
scrublet_result = run_scrublet(adata_list)


##################################################################

# SAVE LAYER OF RAW, FILTERED DATA FOR EACH ANNDATA OBJECT


def save_raw_to_layer(adata_list, layer_name='raw_filtered'):
    for adata in adata_list:
        adata.layers[layer_name] = adata.X.copy()


adata_list = [adata1, adata2, adata3, adata4, adata5]
save_raw_to_layer(adata_list)



##################################################################

# CHECKPOINT A: PRE-SCRUBLET REMOVAL

import copy

# Create a deep copy of each AnnData object to use as a checkpoint
checkpoint_A = tuple([copy.deepcopy(adata) for adata in adata_list])

# Restore adata list from the checkpoint
#adata_list = [copy.deepcopy(adata) for adata in checkpoint_A]

##################################################################

# FILTER DOUBLETS


def filter_doublets(adata_list, doublet_score_threshold=0.5):
    for adata in adata_list:
        # Filter out cells with a doublet score higher than the threshold
        adata._inplace_subset_obs(adata.obs['doublet_scores'] <= doublet_score_threshold)

# Apply the doublet filtering step
adata_list = [adata1, adata2, adata3, adata4, adata5]
filter_doublets(adata_list)


##################################################################

#QC 4

##################################################################

# NORMALIZE, LOGARITHMIZE, SCALE DATA; SAVE LAYERS


def preprocess_adata(adata, layers_dict):
    """
    Preprocess an AnnData object: normalization, logarithmization, and scaling.
    Save each step to a layer.
    """
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

adata_list = [adata1, adata2, adata3, adata4, adata5]
process_adata_list(adata_list)

##################################################################


# QC 4: POST-PROCESSING QC

##################################################################


# CHECKPOINT A - Save h5ad files

#os.chdir("/Users/maureen/Documents/data-experiments/merfish/h5ad-files/checkpoint_A/")

# Replace 'your_filename.h5ad' with the desired file name
#sc.write('adata1_post_scaling.h5ad', adata1)
#sc.write('adata2_post_scaling.h5ad', adata2)
#sc.write('adata3_post_scaling.h5ad', adata3)
#sc.write('adata4_post_scaling.h5ad', adata4)
#sc.write('adata5_post_scaling.h5ad', adata5)




##################################################################

# UMAP EMBEDDINGS


##################################################################


# PCA