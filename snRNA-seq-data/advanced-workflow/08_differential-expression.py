#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 13:16:23 2024

@author: maureen

Notes
- Split Mixed neuronal population
"""

import os
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import diffxpy.api as de
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from functools import reduce
from functools import reduce

plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'Arial'

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '6-tbi-annotated-full.h5ad'))

sc.pl.umap(adata, color = ['cell_type', 'cluster' ,'Slc17a6'])

###############################################################################

# PREPARE ANNDATA

## Convert to dense array
type(adata.X)
adata.X = adata.X.toarray()
type(adata.X)

###############################################################################

# PARAMETERIZED FUNCTION

## size factors function
def calculate_size_factors(subset_data):
    total_counts = subset_data.X.sum(axis=1)
    size_factors = total_counts / np.mean(total_counts)
    return size_factors

## DE function
def run_de_analysis(adata, cell_type, reference_group, comparison_group, min_cells_per_group=40, max_imbalance_ratio=0.9):
    subset_data = adata[(adata.obs['cell_type'] == cell_type) & 
                        (adata.obs['group'].isin([reference_group, comparison_group]))].copy()

    group_counts = subset_data.obs['group'].value_counts()
    if any(group_counts < min_cells_per_group):
        return f"Skipped: {reference_group} vs {comparison_group} (n is too low)"

    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        return f"Skipped: {reference_group} vs {comparison_group} (imbalance)"

    subset_data.X = subset_data.layers['counts'].copy()

    sc.pp.filter_genes(subset_data, min_cells=3)

    if issparse(subset_data.X):
        subset_data.X = subset_data.X.toarray()

    size_factors = calculate_size_factors(subset_data)
    
    subset_data.obs['size_factors'] = size_factors

    test_result = de.test.wald(
        data=subset_data,
        formula_loc="~ 1 + group",
        factor_loc_totest="group",
        size_factors=size_factors
    )
    
    df_result = test_result.summary()
    df_result['cell_type'] = cell_type
    df_result['comparison'] = f"{reference_group}{comparison_group}"
    
    return df_result

###############################################################################

# RUN DE

cell_types = ['Excitatory neuron', 'Inhibitory neuron', 'Mixed neuron', 'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia', 'Choroid-plexus epithelial', 'Fibroblast']

#comparisons = [('A', 'C'), ('C', 'E'), ('A', 'B')]

comparisons = [('A', 'C'), ('C', 'E'), ('A', 'B'), ('A', 'D'), ('C', 'D'), ('D', 'F'), ('E', 'F')]

all_results = []

for cell_type in cell_types:
    for reference_group, comparison_group in comparisons:
        result = run_de_analysis(
            adata=adata,
            cell_type=cell_type,
            reference_group=reference_group,
            comparison_group=comparison_group
        )
        if isinstance(result, pd.DataFrame):
            all_results.append(result)
        else:
            print(f"Skipping {cell_type} {reference_group} vs {comparison_group}: {result}")  

final_results = pd.concat(all_results, ignore_index=True)

###############################################################################

# FILTER RESULTS

def filter_de_results(df, max_log2fc=100):
    return df[(abs(df['log2fc']) <= max_log2fc)]


filtered_results = filter_de_results(final_results)

###############################################################################

# EXPORT DE RESULTS

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/results/de/cell_type"

final_results.to_csv(os.path.join(save_dir, 'wald-test-cell_type.csv'), index = False)

filtered_results.to_csv(os.path.join(save_dir, 'wald-test-cell_type-filtered.csv'), index = False)

###############################################################################

# IMPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/results/de/cell_type"


final_results =  pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type.csv'))

filtered_results =  pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type-filtered.csv'))


save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/results/de/cell_type"

## Export filtered versions for volcano plots
comparison_ids = filtered_results['comparison'].unique()

## Loop through each comparison and save the filtered data
for comp in comparison_ids:
    df_filtered = filtered_results[filtered_results['comparison'] == comp]    
    file_name = f"{comp}-comparison-cell_type-filtered.csv"
    
    # Save the DataFrame as a CSV file
    df_filtered.to_csv(os.path.join(save_dir, file_name), index=False)

###############################################################################

# VISUALIZE

comparisons = ['AC', 'CE', 'AB']

order = ['Astrocyte', 'Excitatory neuron', 'Inhibitory neuron',  'Microglia', 'OPC', 'Oligodendrocyte', 'Fibroblast', 'Choroid-plexus epithelial' ]

## Create a reference DataFrame for cell types
cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

data_list_up = []
data_list_down = []

mean_expression_threshold = 0.1
## Loop over each comparison
for comp in comparisons:
    df_comp = filtered_results[(filtered_results['comparison'] == comp)]
    df_comp = df_comp[df_comp['mean'] > mean_expression_threshold]

    # Upregulated genes
    df_up = df_comp[(df_comp['qval'] < 0.05) & (df_comp['log2fc'] > 0.5)]
    deg_counts_up = df_up.groupby('cell_type').size().reset_index(name=comp)  # Retain original comparison name
    deg_counts_up = pd.merge(cell_type_reference, deg_counts_up, on='cell_type', how='left').fillna(0)
    data_list_up.append(deg_counts_up)

    # Downregulated genes
    df_down = df_comp[(df_comp['qval'] < 0.05) & (df_comp['log2fc'] < -0.5)]
    deg_counts_down = df_down.groupby('cell_type').size().reset_index(name=comp)  # Retain original comparison name
    deg_counts_down = pd.merge(cell_type_reference, deg_counts_down, on='cell_type', how='left').fillna(0)
    data_list_down.append(deg_counts_down)

## Combine up and down dfs
combined_df_up = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list_up)
combined_df_down = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list_down)
combined_df_up.fillna(0, inplace=True)
combined_df_down.fillna(0, inplace=True)

## Calculate vmin and vmax for 
vmin = min(combined_df_up.drop(columns='cell_type').min().min(), combined_df_down.drop(columns='cell_type').min().min())
vmax = max(combined_df_up.drop(columns='cell_type').max().max(), combined_df_down.drop(columns='cell_type').max().max())

# Reorder by cell type
combined_df_up['cell_type'] = pd.Categorical(combined_df_up['cell_type'], categories=order, ordered=True)
combined_df_up = combined_df_up.sort_values('cell_type')

combined_df_down['cell_type'] = pd.Categorical(combined_df_down['cell_type'], categories=order, ordered=True)
combined_df_down = combined_df_down.sort_values('cell_type')


save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/submission_01/figures'

# Plot the heatmap for upregulated genes using the common scale
plt.figure(figsize=(6, 3))
sns.heatmap(combined_df_up.set_index('cell_type'), annot=True, cmap='viridis', fmt="g", vmin=vmin, vmax=vmax)
plt.title('Upregulated Genes')
plt.ylabel(' ')
plt.xlabel('')  # No label for the x-axis
plt.xticks(ticks=range(len(comparisons)), labels=comparisons)  # Set custom x-axis labels
plt.grid(False)
plt.tight_layout()
plt.xticks(ticks=np.arange(len(comparisons)) + 0.5, labels=comparisons, rotation=0)
plt.savefig(os.path.join(save_dir, 'upregulated-genes-summary.pdf'))
plt.show()

# Plot the heatmap for downregulated genes using the common scale
plt.figure(figsize=(6, 3))
sns.heatmap(combined_df_down.set_index('cell_type'), annot=True, cmap='viridis', fmt="g", vmin=vmin, vmax=vmax)
plt.title('Downregulated Genes')
plt.ylabel(' ')
plt.xlabel('')  # No label for the x-axis
plt.xticks(ticks=range(len(comparisons)), labels=comparisons)  # Set custom x-axis labels
plt.grid(False)
plt.tight_layout()
plt.xticks(ticks=np.arange(len(comparisons)) + 0.5, labels=comparisons, rotation=0)
plt.savefig(os.path.join(save_dir, 'downregulated-genes-summary.pdf'))
plt.show()

###############################################################################