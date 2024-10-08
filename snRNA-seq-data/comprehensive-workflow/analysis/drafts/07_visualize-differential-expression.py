#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 16:18:39 2024

@author: maureen
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc
from functools import reduce



plt.rcParams['font.family'] = 'Arial'
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400, fontsize=15)

os.chdir('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/composition')

################################################################################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '3-tbi-annotated-hvg.h5ad'))
adata.X = adata.layers['log1p'].copy()

################################################################################################################################

# PREPARE DATA
save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/composition'

## Exclude artifacts
adata = adata[~adata.obs['cell_type'].isin(['Artifact', 'Unknown'])].copy()

## Update annotations
adata.obs['cell_type'] = adata.obs['cell_type'].replace('Unassigned', 'Unknown')
adata.obs['cell_class'] = adata.obs['cell_class'].replace('Unassigned', 'Unknown')
adata.obs['cluster'] = adata.obs['cluster'].replace('Unassigned', 'Unknown')

################################################################################################################################

# DIFFXPY WALD TEST RESULTS

################################################################################################################################

# WALD TEST - 0.5 LOG2FC CUTOFF

wald_results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de/cell_type/filtered'

# List of comparison identifiers
comparisons = ['AC', 'CE', 'AB']

data_list = []

## Specify order
order = ['Astrocyte', 'Excitatory neuron', 'Fibroblast', 'Inhibitory neuron', 'Microglia', 'OPC', 'Oligodendrocyte', 'Unassigned']

cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

## Loop over each comparison
for comp in comparisons:
    ## Load the filtered results CSV
    file_path = os.path.join(wald_results_dir, f"{comp}-comparison-cell_type-filtered.csv")
    df = pd.read_csv(file_path)
    ## Filter for significant DEGs with a log2FC cutoff
    df_significant = df[(df['qval'] < 0.05) & (abs(df['log2fc']) > 0.5)]
    ## Count number of DEGs by cell type
    deg_counts = df_significant.groupby('cell_type').size().reset_index(name=f'{comp}')
    ## Merge with reference to ensure all cell types are present
    deg_counts = pd.merge(cell_type_reference, deg_counts, on='cell_type', how='left').fillna(0)
    ## Append to  list
    data_list.append(deg_counts)

## Merge all dataframes on 'cell_type'
wald_results_df = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list)
wald_results_df.fillna(0, inplace=True)  # Replace NaN with 0 where there are no DEGs

## Reorder the dataframe according to the desired cell type order
wald_results_df['cell_type'] = pd.Categorical(wald_results_df['cell_type'], categories=order, ordered=True)
wald_results_df = wald_results_df.sort_values('cell_type')

## Rename "Unassigned" to "Unknown" for visualization
wald_results_df['cell_type'] = wald_results_df['cell_type'].replace('Unassigned', 'Unknown')


## Plotting

fig_dir ='/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/de-heatmap'

# Create the heatmap
plt.figure(figsize=(7, 5))
sns.heatmap(wald_results_df.set_index('cell_type'), annot=True, cmap='viridis', fmt="g")
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'wald-deg-heatmap-log2fc-0.5.pdf'))
plt.show()

################################################################################################################################

# WALD TEST - 0.1 log2fc cutoff

wald_results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de/cell_type/filtered'

## List of comparison identifiers
comparisons = ['AC', 'CE', 'AB']

data_list = []

## Specify order
order = ['Astrocyte', 'Excitatory neuron', 'Fibroblast', 'Inhibitory neuron', 'Microglia', 'OPC', 'Oligodendrocyte', 'Unassigned']

cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

## Loop over each comparison
for comp in comparisons:
    ## Load the filtered results CSV
    file_path = os.path.join(wald_results_dir, f"{comp}-comparison-cell_type-filtered.csv")
    df = pd.read_csv(file_path)
    ## Filter for significant DEGs with a log2FC cutoff
    df_significant = df[(df['qval'] < 0.05) & (abs(df['log2fc']) > 0.1)]
    ## Count number of DEGs by cell type
    deg_counts = df_significant.groupby('cell_type').size().reset_index(name=f'{comp}')
    ## Merge with reference to ensure all cell types are present
    deg_counts = pd.merge(cell_type_reference, deg_counts, on='cell_type', how='left').fillna(0)
    ## Append to  list
    data_list.append(deg_counts)

## Merge all dataframes on 'cell_type'
wald_results_df = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list)
wald_results_df.fillna(0, inplace=True)  # Replace NaN with 0 where there are no DEGs

## Reorder the dataframe according to the desired cell type order
wald_results_df['cell_type'] = pd.Categorical(wald_results_df['cell_type'], categories=order, ordered=True)
wald_results_df = wald_results_df.sort_values('cell_type')

## Rename "Unassigned" to "Unknown" for visualization
wald_results_df['cell_type'] = wald_results_df['cell_type'].replace('Unassigned', 'Unknown')


## Plotting

fig_dir ='/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/de-heatmap'

# Create the heatmap
plt.figure(figsize=(7, 5))
sns.heatmap(wald_results_df.set_index('cell_type'), annot=True, cmap='viridis', fmt="g")
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'wald-deg-heatmap-log2fc-0.1.pdf'))
plt.show()


################################################################################################################################

# WALD TEST - 1 log2fc cutoff

wald_results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de/cell_type/filtered'

# List of comparison identifiers
comparisons = ['AC', 'CE', 'AB']

data_list = []

## Specify order
order = ['Astrocyte', 'Excitatory neuron', 'Fibroblast', 'Inhibitory neuron', 'Microglia', 'OPC', 'Oligodendrocyte', 'Unassigned']

cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

## Loop over each comparison
for comp in comparisons:
    ## Load the filtered results CSV
    file_path = os.path.join(wald_results_dir, f"{comp}-comparison-cell_type-filtered.csv")
    df = pd.read_csv(file_path)
    ## Filter for significant DEGs with a log2FC cutoff
    df_significant = df[(df['qval'] < 0.05) & (abs(df['log2fc']) > 1.0)]
    ## Count number of DEGs by cell type
    deg_counts = df_significant.groupby('cell_type').size().reset_index(name=f'{comp}')
    ## Merge with reference to ensure all cell types are present
    deg_counts = pd.merge(cell_type_reference, deg_counts, on='cell_type', how='left').fillna(0)
    ## Append to  list
    data_list.append(deg_counts)

## Merge all dataframes on 'cell_type'
wald_results_df = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list)
wald_results_df.fillna(0, inplace=True)  # Replace NaN with 0 where there are no DEGs

## Reorder the dataframe according to the desired cell type order
wald_results_df['cell_type'] = pd.Categorical(wald_results_df['cell_type'], categories=order, ordered=True)
wald_results_df = wald_results_df.sort_values('cell_type')

## Rename "Unassigned" to "Unknown" for visualization
wald_results_df['cell_type'] = wald_results_df['cell_type'].replace('Unassigned', 'Unknown')


## Plotting

fig_dir ='/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/de-heatmap'

# Create the heatmap
plt.figure(figsize=(7, 5))
sns.heatmap(wald_results_df.set_index('cell_type'), annot=True, cmap='viridis', fmt="g")
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'wald-deg-heatmap-log2fc-1.pdf'))
plt.show()

################################################################################################################################

# WALD TEST - NO log2fc cutoff

wald_results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de/cell_type/filtered'

# List of comparison identifiers
comparisons = ['AC', 'CE', 'AB']

data_list = []

## Specify order
order = ['Astrocyte', 'Excitatory neuron', 'Fibroblast', 'Inhibitory neuron', 'Microglia', 'OPC', 'Oligodendrocyte', 'Unassigned']

cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

## Loop over each comparison
for comp in comparisons:
    ## Load the filtered results CSV
    file_path = os.path.join(wald_results_dir, f"{comp}-comparison-cell_type-filtered.csv")
    df = pd.read_csv(file_path)
    ## Filter for significant DEGs with a log2FC cutoff
    df_significant = df[(df['qval'] < 0.05) & (abs(df['log2fc']) > 0)]
    ## Count number of DEGs by cell type
    deg_counts = df_significant.groupby('cell_type').size().reset_index(name=f'{comp}')
    ## Merge with reference to ensure all cell types are present
    deg_counts = pd.merge(cell_type_reference, deg_counts, on='cell_type', how='left').fillna(0)
    ## Append to  list
    data_list.append(deg_counts)

## Merge all dataframes on 'cell_type'
wald_results_df = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list)
wald_results_df.fillna(0, inplace=True)  # Replace NaN with 0 where there are no DEGs

## Reorder the dataframe according to the desired cell type order
wald_results_df['cell_type'] = pd.Categorical(wald_results_df['cell_type'], categories=order, ordered=True)
wald_results_df = wald_results_df.sort_values('cell_type')

## Rename "Unassigned" to "Unknown" for visualization
wald_results_df['cell_type'] = wald_results_df['cell_type'].replace('Unassigned', 'Unknown')


## Plotting

fig_dir ='/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/de-heatmap'

# Create the heatmap
plt.figure(figsize=(7, 5))
sns.heatmap(wald_results_df.set_index('cell_type'), annot=True, cmap='viridis', fmt="g")
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'wald-deg-heatmap-all-sig.pdf'))
plt.show()

################################################################################################################################

# SCVI

################################################################################################################################

# SCVI RESULTS - 0.5 LOG2FC CUTOFF

results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de-scvi/cell_type'
comparisons = ['AC', 'CE', 'AB'] 

# Collect all unique cell types from all files to ensure none are dropped
all_cell_types = set()
for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    all_cell_types.update(df['Cell_Type'].unique())

# Convert set to list to keep consistent order
all_cell_types = list(all_cell_types)

data_list = []

for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    print(f"Loaded {comp} with shape {df.shape}")

    df_significant = df[(df['is_de_fdr_0.05'] == True) & (abs(df['lfc_median']) > 0.5)]
    print(f"Filtered {comp} significant DEGs with shape {df_significant.shape}")

    deg_counts = df_significant.groupby('Cell_Type').size().reset_index(name=f'{comp}')
    # Initialize a DataFrame for all cell types with zeros
    full_deg_counts = pd.DataFrame({'Cell_Type': all_cell_types})
    full_deg_counts = full_deg_counts.merge(deg_counts, on='Cell_Type', how='left').fillna(0)
    
    data_list.append(full_deg_counts)

# Merge all dataframes on 'Cell_Type'
final_df = reduce(lambda left, right: pd.merge(left, right, on='Cell_Type', how='outer'), data_list)
final_df.fillna(0, inplace=True)

print(final_df)

################################################################################################################################


## Plotting

fig_dir ='/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/scvi/cell_type/updated'


subset_df = final_df[['Cell_Type', 'AC', 'CE', 'AB']]
sns.set_context('talk')

plt.figure(figsize=(7, 5))
sns.heatmap(subset_df.set_index('Cell_Type'), annot=True, cmap='viridis', fmt='g')
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'scvi-deg-heatmap-main-0.5-updated.pdf'))
plt.show()

################################################################################################################################

# SCVI RESULTS - 1.0 LOG2FC CUTOFF

results_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de-scvi/cell_type'
comparisons = ['AC', 'CE', 'AB'] 

# Collect all unique cell types from all files to ensure none are dropped
all_cell_types = set()
for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    all_cell_types.update(df['Cell_Type'].unique())

# Convert set to list to keep consistent order
all_cell_types = list(all_cell_types)

data_list = []

for comp in comparisons:
    file_path = os.path.join(results_dir, f"{comp}-comparison.csv")
    df = pd.read_csv(file_path)
    print(f"Loaded {comp} with shape {df.shape}")

    df_significant = df[(df['is_de_fdr_0.05'] == True) & (abs(df['lfc_median']) > 1)]
    print(f"Filtered {comp} significant DEGs with shape {df_significant.shape}")

    deg_counts = df_significant.groupby('Cell_Type').size().reset_index(name=f'{comp}')
    # Initialize a DataFrame for all cell types with zeros
    full_deg_counts = pd.DataFrame({'Cell_Type': all_cell_types})
    full_deg_counts = full_deg_counts.merge(deg_counts, on='Cell_Type', how='left').fillna(0)
    
    data_list.append(full_deg_counts)

# Merge all dataframes on 'Cell_Type'
final_df = reduce(lambda left, right: pd.merge(left, right, on='Cell_Type', how='outer'), data_list)
final_df.fillna(0, inplace=True)

print(final_df)


## Plotting

fig_dir ='/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/scvi/cell_type/updated'


subset_df = final_df[['Cell_Type', 'AC', 'CE', 'AB']]
sns.set_context('talk')

plt.figure(figsize=(7, 5))
sns.heatmap(subset_df.set_index('Cell_Type'), annot=True, cmap='viridis', fmt='g')
plt.title('')
plt.ylabel('')
plt.xlabel('')
plt.grid(False)
plt.tight_layout()
plt.savefig(os.path.join(fig_dir, 'scvi-deg-heatmap-main-1.0-updated.pdf'))
plt.show()

################################################################################################################################