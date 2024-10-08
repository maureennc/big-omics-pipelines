#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 14:06:04 2024

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
from functools import reduce

import random
import torch
import scvi


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
plt.rcParams['figure.dpi'] = 300

## Scanpy
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

################################################################################################################################

# IMPORT DATA

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '3-tbi-annotated-hvg.h5ad'))
adata.X = adata.layers['log1p'].copy()
adata.obs['cell_type'] = adata.obs['cell_type'].replace('Unassigned', 'Unknown') # Updated nomenclature 6/24/24 for data viz



## scVI model
scvi_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/rivanna/scvi"
model_dir = os.path.join(scvi_dir, 'model_H')
print(model_dir)

model_H = scvi.model.SCVI.load(model_dir, adata=adata)
model_H

################################################################################################################################

# FUNCTION FOR DE RESULTS

def run_scvi_de(cell_type, group1, group2, model):
    idx1 = (adata.obs['cell_type'] == cell_type) & (adata.obs['group'] == group1)
    idx2 = (adata.obs['cell_type'] == cell_type) & (adata.obs['group'] == group2)
    result = model.differential_expression(idx1=idx1, idx2=idx2)
    result['Cell_Type'] = cell_type  # Add a column for the cell type
    result.reset_index(inplace=True)  # Reset index to bring gene names into a column
    return result

## Setup directory
save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de-scvi/cell_type'

cell_types = ['Microglia', 'Astrocyte', 'Excitatory neuron', 'Inhibitory neuron', 'Oligodendrocyte', 'OPC', 'Fibroblast', 'Unknown']
group_comparisons = [('A', 'C'), ('C', 'E'), ('A', 'B'), ('A', 'D'), ('C', 'D'), ('D', 'F')]

## Collect results
for group1, group2 in group_comparisons:
    comparison_results = []
    for cell_type in cell_types:
        de_results = run_scvi_de(cell_type, group1, group2, model_H)
        comparison_results.append(de_results)
    ## Concatenate results for all cell types for the current comparison
    combined_results = pd.concat(comparison_results, ignore_index=True)
    file_name = f'{group1}{group2}-comparison.csv'
    full_path = os.path.join(save_dir, file_name)
    combined_results.to_csv(full_path, index=False)  # Ensure index=False is set if the gene names are already a column
    print(f"Saved: {full_path}")



################################################################################################################################

# ALL CELL TYPES, ALL GROUPS

results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de-scvi/cell_type'
comparisons = ['AC', 'CE', 'AB', 'AD', 'CD', 'DF'] 

## Collect all unique cell types from all files to ensure none are dropped
all_cell_types = set()
for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    all_cell_types.update(df['Cell_Type'].unique())

## Convert set to list to keep consistent order
all_cell_types = list(all_cell_types)

data_list = []

for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    print(f"Loaded {comp} with shape {df.shape}")

    df_significant = df[df['is_de_fdr_0.05'] & (abs(df['lfc_median']) > 0.5)] # LOG2FC THRESHOLD
    print(f"Filtered {comp} significant DEGs with shape {df_significant.shape}")

    deg_counts = df_significant.groupby('Cell_Type').size().reset_index(name=f'{comp}')
    # Initialize a DataFrame for all cell types with zeros
    full_deg_counts = pd.DataFrame({'Cell_Type': all_cell_types})
    full_deg_counts = full_deg_counts.merge(deg_counts, on='Cell_Type', how='left').fillna(0)
    
    data_list.append(full_deg_counts)

## Merge all dataframes on 'Cell_Type'
final_df = reduce(lambda left, right: pd.merge(left, right, on='Cell_Type', how='outer'), data_list)
final_df.fillna(0, inplace=True)

print(final_df)


fig_dir ='/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/scvi/cell_type'
sns.set_context('talk')

plt.figure(figsize=(7, 5))
sns.heatmap(final_df.set_index('Cell_Type'), annot=True, cmap='viridis', fmt='g')
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'scvi-deg-heatmap-all-0.5.pdf'))
plt.show()


################################################################################################################################

# PRIORITY SUBSET

subset_df = final_df[['Cell_Type', 'AC', 'CE', 'AB']]
sns.set_context('talk')

plt.figure(figsize=(7, 5))
sns.heatmap(subset_df.set_index('Cell_Type'), annot=True, cmap='viridis', fmt='g')
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'scvi-deg-heatmap-main-0.5.pdf'))
plt.show()

################################################################################################################################
