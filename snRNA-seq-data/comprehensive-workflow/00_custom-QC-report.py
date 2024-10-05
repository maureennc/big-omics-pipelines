#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 13 11:37:38 2024

@author: maureen
"""


import os
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
import pandas as pd
import scanpy as sc

# environment = 'sc-pp'

###############################################################################

# RAW DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/azenta-files/01_analysis/cellranger_count"

sample_A = sc.read_10x_mtx(os.path.join(data_dir, "1/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_B = sc.read_10x_mtx(os.path.join(data_dir, "3/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_C = sc.read_10x_mtx(os.path.join(data_dir, "5/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_D = sc.read_10x_mtx(os.path.join(data_dir, "6/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_E = sc.read_10x_mtx(os.path.join(data_dir, "7/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_F = sc.read_10x_mtx(os.path.join(data_dir, "8/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)

adata_list = [sample_A, sample_B, sample_C, sample_D, sample_E, sample_F]

for sample in adata_list:
    print(sample)
    
###############################################################################

# CONCATENATION

## Perform concatenation
adata_concat = sample_A.concatenate(sample_B, sample_C, sample_D, sample_E, sample_F, batch_key='group', batch_categories=['A', 'B', 'C', 'D', 'E', 'F'] )

# QC
adata_concat.var['mt'] = adata_concat.var_names.str.startswith('mt-')
adata_concat.var['ribosomal'] = adata_concat.var_names.str.match('^(Rpl|Rps)\\d+')

adata_concat.obs['total_reads_per_cell'] = adata_concat.X.sum(axis=1)
sc.pp.calculate_qc_metrics(adata_concat, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(adata_concat, keys = 'total_reads_per_cell', groupby = 'group')
sc.pl.violin(adata_concat, keys = 'total_counts', groupby = 'group', ylabel = '# Transcripts per cell')

## Calculate mean total_counts
mean_total_counts_by_sample = adata_concat.obs.groupby('group')['total_counts'].mean()
mean_total_counts_by_sample

## Calculate mean genes per cell
mean_unique_genes = adata_concat.obs.groupby('group')['n_genes_by_counts'].mean()
mean_unique_genes
    
###############################################################################

# QC-PASSED DATA

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(save_dir, '1-tbi-seq-hvg.h5ad'))
adata.X = adata.layers['counts'].copy()

sc.pl.violin(adata, keys = 'total_counts', groupby = 'group', ylabel = "# Transcripts per cell")

## Calculate median total_counts
mean_total_counts_by_sample = adata.obs.groupby('group')['total_counts'].mean()
mean_total_counts_by_sample

## Calculate mean genes per cell
mean_unique_genes = adata.obs.groupby('group')['n_genes_by_counts'].mean()
mean_unique_genes
    
###############################################################################

# DOUBLET PASSED DATA

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

bdata = sc.read_h5ad(os.path.join(save_dir, '3-tbi-annotated-hvg.h5ad'))

bdata.obs.cell_type.value_counts()

sc.pl.violin(bdata, keys = 'total_counts', groupby = 'group', ylabel = "# Transcripts per cell")

## Calculate median total_counts
median_total_counts_by_sample = bdata.obs.groupby('group')['total_counts'].median()
median_total_counts_by_sample

mean_total_counts_by_sample = bdata.obs.groupby('group')['total_counts'].mean()
mean_total_counts_by_sample

## Calculate mean genes per cell
mean_unique_genes = bdata.obs.groupby('group')['n_genes_by_counts'].mean()
mean_unique_genes

###############################################################################

# CREATE SUMMARY DF

data = {
    'group': ['A', 'B', 'C', 'D', 'E', 'F'],
    '# cells cellranger_filtered': [6765, 10109, 7926, 9132, 7233, 5945],
    '# cells pass initial QC': [5144, 7196, 5765, 6574, 5304, 4325],
    '# after doublet detection': [4355, 6063, 5011, 5475, 4304, 3846]
}

# Convert  dictionary into a pandas df
df = pd.DataFrame(data)

# Display the DataFrame
print(df)

## Add percent pass for each
## Add percent columns to the DataFrame
df['% pass initial QC'] = (df['# cells pass initial QC'] / df['# cells cellranger_filtered']) * 100
df['% pass doublet detection'] = (df['# after doublet detection'] / df['# cells pass initial QC']) * 100

df['% high-quality'] = (df['# after doublet detection'] / df['# cells cellranger_filtered']) * 100

print(df)

###############################################################################

# SUMMARIZE DOUBLET RESULTS FROM SOLO

solo_dir ="/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/solo"

solo = pd.read_csv(os.path.join(solo_dir, 'solo-results-dr2.csv'))

counts = solo.prediction.value_counts()

total_cells = len(solo)

# Calculate the proportion of singlets and doublets
singlet_proportion = counts['singlet'] / total_cells
doublet_proportion = counts['doublet'] / total_cells

print(f"Estimated proportion of singlets: {singlet_proportion:.2%}")
print(f"Estimated proportion of doublets: {doublet_proportion:.2%}")

###############################################################################