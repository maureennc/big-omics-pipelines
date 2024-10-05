#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 23 12:40:56 2024

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
import scipy.stats
from scipy.sparse import csr_matrix
import diffxpy.api as de
import random


################################################################################################################################

# SETTINGS

## Random seed
random.seed(0)
np.random.seed(0)

## Matplotlib
%matplotlib inline
plt.rcParams['font.family'] = 'Arial'

## Scanpy
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, 'mg-model_B-cleaned-annotated.h5ad'))


## Save directory
save_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/mg-sting-ko/results/de"


################################################################################################################################

# PREPARE DATA

adata.X = adata.X.toarray()

## Restore the raw counts for updated normalization
adata.X = adata.layers['counts'].copy()

## Apply median-based normalization
sc.pp.normalize_total(adata)
adata.layers['normalized'] = adata.X.copy()

## Apply log1p transformation
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()

## Freeze log1p counts in raw slot
adata.raw = adata.copy()

################################################################################################################################

# PERFORM DE AT CELL TYPE LEVEL

results_list = []
cell_types = adata.obs['cell_type'].unique()

## Loop through each  cell type and prepare the data
for cell_type in cell_types:
    cell_type_data = adata[adata.obs['cell_type'] == cell_type].copy()
    groups = cell_type_data.obs['group'].value_counts()
    
    # Filter groups having at least 50 cells
    groups = groups[groups >= 40].index.tolist()
    
    # Ensure there are sufficient cells and groups meet the criteria
    if cell_type_data.shape[0] > 100 and len(groups) > 1:
        # Check for balance: no group should be less than 10% of the cell type size
        balanced = all(cell_type_data.obs['group'].value_counts(normalize=True) > 0.10)
        
        if balanced:
            cell_type_data.X = cell_type_data.X.toarray()  # Convert sparse matrix to dense array if necessary
            try:
                # Perform DE
                test_result = de.test.wald(
                    data=cell_type_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cell_type'] = cell_type
                results_list.append(df_result)
                print(f"Processed cell type: {cell_type}")
            except Exception as e:
                print(f"Error processing cell type {cell_type}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cell_type} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cell_type} due to insufficient data or group criteria not met.")

## Concatenate
if results_list:
    results_cell_type = pd.concat(results_list, ignore_index=True)
    print(results_cell_type)
else:
    print("No results to display.")
    
 
results_cell_type = pd.concat(results_list, ignore_index=True)

filtered_de_cell_type = results_cell_type[(abs(results_cell_type['log2fc']) <= 100) & (results_cell_type['mean'] >= 0.05)]
filtered_de_cell_type['log2fc'] = -filtered_de_cell_type['log2fc']
#filtered_de_cell_type.to_csv(os.path.join(save_dir, 'mg-sting-wt-ko-DE-cell_type-median-norm.csv'), index = False)

################################################################################################################################

# VERIFY REFRENCE VS. TESTING GROUP

specific_cell_type = 'Astrocyte'  # Change this to cell type of interest

## Filter the DataFrame for the specific cell type
de_df_specific = filtered_de_cell_type[filtered_de_cell_type['cell_type'] == specific_cell_type]

## Sort by log2fc to find the most upregulated genes, descending
de_df_sorted = de_df_specific.sort_values('log2fc', ascending=False).reset_index(drop=True)

if not de_df_sorted.empty:
    most_upregulated = de_df_sorted.iloc[0].gene
    print(f"The most upregulated gene in {specific_cell_type} is: {most_upregulated}")
    
    # Get the index of the most upregulated gene in the AnnData var_names
    i = np.where(adata.var_names == most_upregulated)[0][0]
    
    # Extract expression data for the most upregulated gene from two groups within the specified cell type
    filtered_adata = adata[(adata.obs['group'].isin(['WT', 'Sting-KO'])) & (adata.obs['cell_type'] == specific_cell_type)]
    a = filtered_adata[filtered_adata.obs['group'] == 'WT', i].X
    b = filtered_adata[filtered_adata.obs['group'] == 'Sting-KO', i].X
    
    # Calculate the mean expression in each group
    a_mean = np.mean(a.toarray()) if hasattr(a, "toarray") else np.mean(a)
    b_mean = np.mean(b.toarray()) if hasattr(b, "toarray") else np.mean(b)
    
    print(f"{most_upregulated} expression in {specific_cell_type}:")
    print(f"WT: {a_mean}")
    print(f"Sting-KO: {b_mean}")
else:
    print(f"No differential expression results available for {specific_cell_type}")

# Doesn't work for endotheilal; gene is Ifitm1

################################################################################################################################

# PERFORM DE AT CELL CLUSTER LEVEL

results_list = []
cluster = adata.obs['cluster'].unique()

# Loop through each cell type and prepare the data
for cluster in cluster:
    cluster_data = adata[adata.obs['cluster'] == cluster].copy()
    groups = cluster_data.obs['group'].value_counts()
    
    groups = groups[groups >= 40].index.tolist()
    
    if cluster_data.shape[0] > 100 and len(groups) > 1:
        balanced = all(cluster_data.obs['group'].value_counts(normalize=True) > 0.10)
        if balanced:
            cluster_data.X = cluster_data.X.toarray()
            try:
                test_result = de.test.wald(
                    data=cluster_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cluster'] = cluster
                results_list.append(df_result)
                print(f"Processed cell type: {cluster}")
            except Exception as e:
                print(f"Error processing cell type {cluster}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cluster} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cluster} due to insufficient data or group criteria not met.")

# Concatenate
if results_list:
    results_cluster = pd.concat(results_list, ignore_index=True)
    print(results_cluster)
else:
    print("No results to display.")
    
 
results_cluster = pd.concat(results_list, ignore_index=True)

filtered_de_cluster = results_cluster[(abs(results_cluster['log2fc']) <= 100) & (results_cluster['mean'] >= 0.05)]
filtered_de_cluster['log2fc'] = -filtered_de_cluster['log2fc']
filtered_de_cluster.to_csv(os.path.join(save_dir, 'mg-sting-wt-ko-DE-cluster-median-norm.csv'), index = False)

################################################################################################################################