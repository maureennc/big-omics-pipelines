#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:14:07 2024

@author: maureen
"""

# CITE-SEQ PRE-PROCESSING

## env = 'sc-pp'

################################################################################################################################

# IMPORT PACKAGES AND SETTINGS

import matplotlib
import matplotlib.pyplot as plt
import scanpy as sc
import os
import anndata as ad
import seaborn as sns
import scrublet as scr
import pandas as pd
from scipy.sparse import csr_matrix


%matplotlib qt5
plt.rcParams['font.family'] = 'Arial'

data_dir = "/Users/maureen/Documents/experiments/cite-seq/data/cellranger-output/output-files-velocyto"

################################################################################################################################

# IMPORT DATA

naive = sc.read_10x_mtx(os.path.join(data_dir, "naive/filtered_feature_bc_matrix"), gex_only = False)
INF1 = sc.read_10x_mtx(os.path.join(data_dir, "INF1/filtered_feature_bc_matrix"), gex_only = False)
INF2 = sc.read_10x_mtx(os.path.join(data_dir, "INF2/filtered_feature_bc_matrix"), gex_only = False)

################################################################################################################################

# SEPARATE GENE EXPRESSION AND PROTEIN DATA

## Add protein feature types to adata.obsm
protein_indices = naive.var[naive.var['feature_types'] == 'Antibody Capture'].index
protein_naive = naive[:, protein_indices].X
naive.obsm['protein_expression'] = protein_naive

protein_indices = INF1.var[INF1.var['feature_types'] == 'Antibody Capture'].index
protein_INF1 = INF1[:, protein_indices].X
INF1.obsm['protein_expression'] = protein_INF1

protein_indices = INF2.var[INF2.var['feature_types'] == 'Antibody Capture'].index
protein_INF2 = INF2[:, protein_indices].X
INF2.obsm['protein_expression'] = protein_INF2


## Examine protein data
print(naive.obsm['protein_expression'].shape)
print(INF1.obsm['protein_expression'].shape)
print(INF2.obsm['protein_expression'].shape)


## Slice adata.X to contain only GEX data
naive = naive[:, naive.var["feature_types"] == "Gene Expression"].copy()
INF1 = INF1[:, INF1.var["feature_types"] == "Gene Expression"].copy()
INF2 = INF2[:, INF2.var["feature_types"] == "Gene Expression"].copy()


################################################################################################################################

## CALCULATE MITOCHONDRIAL AND RIBOSOMAL GENES

adata_list = [naive, INF1, INF2]

for i, adata in enumerate(adata_list):
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    adata.var['ribosomal'] = adata.var_names.str.match('^(Rpl|Rps)\\d+')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)
    
################################################################################################################################

# VISUALIZE QC DATA

adata_list = [naive, INF1, INF2]

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

naive, INF1, INF2 = adata_list

## Mitochondrial percentage filter
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before mitochondrial filter threshold: {adata.shape}")
    adata_list[i] = adata[adata.obs.pct_counts_mt < 5, :].copy()
    print(f"Dataset {i} shape after mitochondrial filter threshold: {adata.shape}")

naive, INF1, INF2 = adata_list


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
    
    
################################################################################################################################

# VISUALIZE QC DATA

adata_list = [naive, INF1, INF2]

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])
    

################################################################################################################################

# DOUBLET DETECTION

def run_scrublet(adata_list):

    scrublet_rows = []

    for i, adata in enumerate(adata_list):
        ## Set up Scrublet
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.10, random_state = 0)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                                  min_cells=3, 
                                                                  min_gene_variability_pctl=85, 
                                                                  n_prin_comps=50)
        adata.obs['doublet_scores'] = doublet_scores
        adata.obs['predicted_doublets'] = predicted_doublets

        ## Plot Histograms
        plt.figure(figsize=(10, 6))
        sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
        sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
        plt.title(f'Scrublet Doublet Score Distribution for Sample {i+1}')
        plt.xlabel('Doublet Score')
        plt.ylabel('Density')
        plt.legend()
        plt.grid(False)
        plt.show()

        ## Extract barcodes
        cell_barcodes = adata.obs.index

        ## Store Scrublet data with barcodes for each sample in list
        for barcode, obs_score, sim_score, pred_doublet in zip(cell_barcodes, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets):
            scrublet_rows.append({'Sample_Index': i+1, 
                                  'Cell_Barcode': barcode,
                                  'Observed_Score': obs_score, 
                                  'Simulated_Score': sim_score, 
                                  'Predicted_Doublet': pred_doublet})

    ## Create df from list of rows
    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df

## Run Scrublet
scrublet_result = run_scrublet(adata_list)

################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/3/h5ad"

naive.X = csr_matrix(naive.X)
INF1.X = csr_matrix(INF1.X)
INF2.X = csr_matrix(INF2.X)

naive.write_h5ad(os.path.join(save_dir, '1-naive-qc-filtered.h5ad'))
INF1.write_h5ad(os.path.join(save_dir, '1-INF1-qc-filtered.h5ad'))
INF2.write_h5ad(os.path.join(save_dir, '1-INF2-qc-filtered.h5ad'))

################################################################################################################################


