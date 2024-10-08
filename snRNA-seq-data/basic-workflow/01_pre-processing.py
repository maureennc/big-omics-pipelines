#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 13:45:21 2024

@author: maureen
"""

import os
from pathlib import Path
import re
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
plt.rcParams['font.family'] = 'Arial'

import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
#sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

from scipy.sparse import csr_matrix
import scipy
import scrublet as scr

# environment = 'sc-pp'

################################################################################################################################

# IMPORT DATA

## WT vs. Microglia-cre KO
mg_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/2_5xSting_MG-KO/"

mg_wt = sc.read_10x_mtx(os.path.join(mg_dir, "Sample3_5xSting_MGWT/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
mg_ko = sc.read_10x_mtx(os.path.join(mg_dir, "Sample4_5xSting_MGKO/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)

## Check for Sting1
print("Gene present in adata.var:", 'Sting1' in mg_wt.var.index)
print("Gene present in adata.var:", 'Sting1' in mg_ko.var.index)

print("Gene present in adata.var:", 'Tmem173' in mg_wt.var.index)
print("Gene present in adata.var:", 'Tmem173' in mg_ko.var.index)

################################################################################################################################

## ANNOTATE MITOCHONDRIAL AND RIBOSOMAL GENES

adata_list = [mg_wt, mg_ko]

for i, adata in enumerate(adata_list):
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribosomal'] = adata.var_names.str.match('^(Rpl|Rps)\\d+')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)


################################################################################################################################

# VISUALIZE QC DATA

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])
    

################################################################################################################################

# FILTER CELLS BASED ON QC THRESHOLDS

## Total counts
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['total_counts'] > 1000, :].copy()
    print(f"Dataset {i} shape after total_counts floor threshold: {adata_list[i].shape}")


### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts ceiling threshold: {adata.shape}")
    total_counts_ceiling = adata.obs['total_counts'].quantile(0.95)
    adata_list[i] = adata[adata.obs['total_counts'] < total_counts_ceiling, :]
    print(f"Dataset {i} shape after total_counts ceiling threshold: {adata_list[i].shape}")

## Number of genes per cell
### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts floor threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] > 500, :].copy()
    print(f"Dataset {i} shape after n_genes_by_counts floor threshold: {adata_list[i].shape}")

### Ceiling threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before n_genes_by_counts ceiling threshold: {adata.shape}")
    ngenes_ceiling = adata.obs['n_genes_by_counts'].quantile(0.97)
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] < ngenes_ceiling, :]
    print(f"Dataset {i} shape after n_genes_by_counts ceiling threshold: {adata_list[i].shape}")

## Ribosomal counts
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before pct_counts_ribosomal ceiling threshold: {adata.shape}")
    ribo_ceiling = adata.obs['pct_counts_ribosomal'].quantile(0.99)
    adata_list[i] = adata[adata.obs['pct_counts_ribosomal'] < ribo_ceiling, :]
    print(f"Dataset {i} shape after pct_counts_ribosomal ceiling threshold: {adata_list[i].shape}")

mg_wt, mg_ko = adata_list

## Mitochondrial percentage filter
for i, adata in enumerate(adata_list):
    adata_list[i] = adata[adata.obs.pct_counts_mt < 5, :].copy()

mg_wt, mg_ko = adata_list

################################################################################################################################

# RECALCULATE QC METRICS

for i, adata in enumerate(adata_list):
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

## Print QC summary data
qc_metrics = ['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal', 'n_genes_by_counts']

for i, adata in enumerate(adata_list, start=1):
    print(f"Dataset {i}:")
    for metric in qc_metrics:
        mean_value = adata.obs[metric].mean()
        print(f"Mean {metric}: {mean_value}")
    print("-" * 30) 

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])


################################################################################################################################

# DOUBLET DETECTION

def run_scrublet(adata_list):

    scrublet_rows = []

    for i, adata in enumerate(adata_list):
        # Run Scrublet
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.10, random_state = 0)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=50)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        # Plot Histograms
        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for Sample {i+1}')
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

# Run Scrublet
scrublet_result = run_scrublet(adata_list)


################################################################################################################################

# FILTER DOUBLETS

def filter_doublets(adata_list):
    filtered_list = []
    for i, adata in enumerate(adata_list):
        print(f"Dataset {i + 1} shape before doublet filter: {adata.shape}")
        filtered_adata = adata[adata.obs['predicted_doublets'] == False].copy()
        filtered_list.append(filtered_adata)
    return filtered_list

## Perform filtering
adata_list = filter_doublets(adata_list)

## Update list reference
mg_wt, mg_ko = adata_list

for i, adata in enumerate(adata_list):
    print(f"Dataset {i + 1} shape after doublet filter: {adata.shape}")


################################################################################################################################

# CONCATENATION

## Annotate groups and methods
mg_wt.obs['group'] = 'WT'
mg_ko.obs['group'] = 'Sting-KO'


## Perform concatenation
adata = mg_wt.concatenate(mg_ko, batch_key='dataset', batch_categories=['mg_wt', 'mg_ko'] )
adata.obs

################################################################################################################################

# FINISH PRE-PROCESSING

## Filter genes appearing in few cells in concatenated dataset
sc.pp.filter_genes(adata, min_cells=3)

## Normalization and Log1p transformation after concatenation
adata.layers['counts'] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum = 1e4)
adata.layers['normalized'] = adata.X.copy()
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()
adata.raw = adata

## Create separate object for full gene set
bdata = adata.copy()

################################################################################################################################

# PRELIMINARY SELECTION OF GENES OF INTEREST AND HVGs

## Import genes of interest list
file_path = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/gene-lists/genes-of-interest.csv"
gene_list = pd.read_csv(file_path)

## Initialize the 'gene_of_interest' column to False
adata.var['gene_of_interest'] = False

# Expand gene list with wildcards and annotate adata.var 'gene_of_interest' column
expanded_rows = []
for _, row in gene_list.iterrows():
    if pd.isna(row['expanded']):  # If there's no expansion, include the gene directly
        if row['gene'] in adata.var_names:
            adata.var.loc[row['gene'], 'gene_of_interest'] = True
    else:
        expansions = row['expanded'].split(',')
        for exp in expansions:
            expanded_gene = f"{row['gene'].rstrip('*')}{exp}"
            if expanded_gene in adata.var_names:
                adata.var.loc[expanded_gene, 'gene_of_interest'] = True

## Select Highly Variable Genes (HVGs)
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=False, layer='counts', flavor="seurat_v3")
adata.var.loc[adata.var['highly_variable'], 'gene_of_interest'] = True

# Subset adata to include only genes of interest
adata_subset = adata[:, adata.var['gene_of_interest']].copy()

## Verification
print(adata_subset.var[['highly_variable', 'gene_of_interest']].value_counts())

## Scale data
sc.pp.scale(adata_subset, max_value=10)
adata_subset.layers['scaled'] = adata_subset.X.copy()

sc.pp.neighbors(adata_subset, random_state = 0)
sc.tl.umap(adata_subset)
sc.tl.leiden(adata_subset, resolution=0.5)
sc.pl.umap(adata_subset, color = ['leiden'])


################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/h5ad"

adata_subset.X = csr_matrix(adata_subset.X)
bdata.X = csr_matrix(bdata.X)

adata_subset.write_h5ad(os.path.join(save_dir, 'mg-sting-ko-pp-concat-hvg.h5ad'))
bdata.write_h5ad(os.path.join(save_dir, 'mg-sting-ko-pp-concat-full.h5ad'))

################################################################################################################################

