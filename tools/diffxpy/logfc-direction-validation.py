#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 09:27:00 2024

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

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/nestin-sting-ko/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, 'nestin-model_B-cleaned-annotated.h5ad'))
adata.X = adata.layers['log1p'].copy()
adata.obs['group'] = adata.obs['group'].replace('Nestin-KO', 'Sting-KO')

de_df_cell_type = pd.read_csv('/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/nestin-sting-ko/results/spreadsheets/de/nestin-sting-wt-ko-DE-cell_type.csv')

################################################################################################################################

# VERIFY REFERENCE VS. TESTING GROUP


specific_cell_type = 'Excitatory neuron'  # Change this to cell type of interest

# Filter the DataFrame for the specific cell type
de_df_specific = de_df_cell_type[de_df_cell_type['cell_type'] == specific_cell_type]

# Sort by 'log2fc' to find the most upregulated genes, descending order
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
    # Assuming your data might be stored in a sparse matrix format, hence using `.toarray()`
    a_mean = np.mean(a.toarray()) if hasattr(a, "toarray") else np.mean(a)
    b_mean = np.mean(b.toarray()) if hasattr(b, "toarray") else np.mean(b)
    
    print(f"{most_upregulated} expression in {specific_cell_type}:")
    print(f"WT: {a_mean}")
    print(f"Sting-KO: {b_mean}")
else:
    print(f"No differential expression results available for {specific_cell_type}")
    
################################################################################################################################

