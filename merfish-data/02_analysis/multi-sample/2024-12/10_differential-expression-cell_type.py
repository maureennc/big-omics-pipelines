#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 10:13:02 2024

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
import squidpy as sq
import random

###############################################################################

# SETTINGS


## Matplotlib
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 500

###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))

###############################################################################

# DE

## Reference is naive
## Group is infected

import pandas as pd
import scanpy as sc

def pairwise_comparison(adata, reference_group='Naive', comparison_group='Infected', min_cells=3):
    # Initialize an empty list to store results
    results = []

    # Loop through all unique cell types in the dataset
    for cell_type in adata.obs['cell_type'].unique():
        
        # Subset the data for the specific cell type
        adata_subset = adata[adata.obs['cell_type'] == cell_type].copy()
        
        # Filter genes to ensure valid genes are used in the comparison
        sc.pp.filter_genes(adata_subset, min_cells=min_cells)
        
        # Check if both groups (Naive and Infected) exist and have enough cells
        group_counts = adata_subset.obs['condition'].value_counts()
        if (reference_group in group_counts and 
            comparison_group in group_counts and 
            group_counts[reference_group] > 40 and 
            group_counts[comparison_group] > 40):
            
            print(f"Comparing {reference_group} vs {comparison_group} for cell type: {cell_type}")
            
            # Run Wilcoxon rank-sum test
            sc.tl.rank_genes_groups(
                adata_subset, 
                groupby='condition', 
                reference=reference_group, 
                groups=[comparison_group], 
                method='wilcoxon', 
                layer='log1p', 
                use_raw=False, 
                pts=True
            )
            
            # Extract DE results
            de_results = sc.get.rank_genes_groups_df(adata_subset, group=comparison_group)
            
            # Calculate pct_nz (percentage of non-zero expression) for both groups
            comparison_mask = adata_subset.obs['condition'] == comparison_group
            reference_mask = adata_subset.obs['condition'] == reference_group
            
            comparison_nonzero = (adata_subset[comparison_mask].X > 0).mean(axis=0).A1  # Sparse matrix handling
            reference_nonzero = (adata_subset[reference_mask].X > 0).mean(axis=0).A1

            # Calculate the mean percentage of non-zero expression
            pct_nz_mean = (comparison_nonzero + reference_nonzero) / 2
            
            # Add calculated values to the DE results DataFrame
            de_results['pct_nz_group'] = comparison_nonzero
            de_results['pct_nz_reference'] = reference_nonzero
            de_results['pct_nz_mean'] = pct_nz_mean
            
            # Add metadata columns
            de_results['cell_type'] = cell_type
            de_results['reference_group'] = reference_group
            de_results['comparison_group'] = comparison_group
            de_results['comparison'] = f'{reference_group}_vs_{comparison_group}'
            
            # Store the results
            results.append(de_results)
        
        else:
            print(f"Skipping {cell_type} due to insufficient cells in {reference_group} or {comparison_group}.")
    
    # Combine all results into a single DataFrame
    combined_df = pd.concat(results, ignore_index=True)
    
    return combined_df

# Run the comparison for Naive vs. Infected conditions
naive_vs_infected_df = pairwise_comparison(adata)

# View the top results (optional)
naive_vs_infected_df.head()

###############################################################################

# HEATMAP

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Ensure the DataFrame is ordered according to the specified cell type order
order = [
    'Astrocyte', 'Choroid plexus', 'Endothelial cell', 'Excitatory neuron', 
    'Inhibitory neuron', 'Microglia', 'OPC', 'Oligodendrocyte', 'Pericyte'
]

cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

# Define thresholds for filtering significant DEGs
logfc_threshold = 0.5  # Use logfoldchanges directly
qval_threshold = 0.01
mean_expression_threshold = 0.1  # Adjust this value as needed

# Filter the results DataFrame for significant DEGs with logfoldchanges and mean expression threshold
filtered_results = naive_vs_infected_df[
    (abs(naive_vs_infected_df['logfoldchanges']) > logfc_threshold) &
    (naive_vs_infected_df['pvals_adj'] < qval_threshold) &
    (naive_vs_infected_df['pct_nz_mean'] > mean_expression_threshold)  # New threshold
]

# Initialize DataFrames to store upregulated and downregulated DEG counts
deg_counts_all = cell_type_reference.copy()
deg_counts_all['Downregulated'] = 0  # Start with zero downregulated counts
deg_counts_all['Upregulated'] = 0  # Start with zero upregulated counts

# Loop over each unique cell type and count DEGs
for cell_type in filtered_results['cell_type'].unique():
    # Subset the results for the specific cell type
    df_cell_type = filtered_results[filtered_results['cell_type'] == cell_type]

    # Count upregulated and downregulated DEGs
    upregulated_count = df_cell_type[df_cell_type['logfoldchanges'] > 0].shape[0]
    downregulated_count = df_cell_type[df_cell_type['logfoldchanges'] < 0].shape[0]

    # Update the counts in the summary DataFrame
    deg_counts_all.loc[deg_counts_all['cell_type'] == cell_type, 'Upregulated'] = upregulated_count
    deg_counts_all.loc[deg_counts_all['cell_type'] == cell_type, 'Downregulated'] = downregulated_count

# Reorder the DataFrame so 'Upregulated' is the first row and 'Downregulated' the second
deg_counts_all = deg_counts_all[['cell_type', 'Upregulated', 'Downregulated']]
deg_counts_all = deg_counts_all.set_index('cell_type').T  # Transpose the DataFrame

# Rename the index labels to "UP" and "DOWN"
deg_counts_all.index = ['UP', 'DOWN']

# Plot the horizontal heatmap
plt.figure(figsize=(5, 2))  # Adjust the figure size to be wider
sns.heatmap(
    deg_counts_all, 
    annot=True, 
    cmap='rocket', 
    fmt="g", 
    cbar_kws={'label': '# DEGs'}
)

plt.title('')
plt.ylabel('') 
plt.xlabel('') 
plt.xticks(rotation=45, ha='right')
plt.grid(False)
plt.tight_layout()
plt.show()
