#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 19 13:57:34 2024

@author: maureen

Virtual environment: 'sc-seq'

This script performs the following steps, in order:
    1. Import cellranger output files from custom genome build (for downstream trajectory inference analysis)
    2. Separate and prepare protein and gene expression data
    2. Filtering and QC
    3. Doublet detection and removal
    4. Concatenation
    5. Normalization
    6. Highly-variable gene selection for dimensionality reduction
    7. Exports processed AnnData as 'cite-seq-pp-concat-hvg.h5ad'
    8. Exports processed AnnData as 'cite-seq-pp-concat-full.h5ad'

"""

# CITE-SEQ PRE-PROCESSING

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

###############################################################################

# IMPORT DATA

naive = sc.read_10x_mtx(os.path.join(data_dir, "naive/filtered_feature_bc_matrix"), gex_only = False)
INF1 = sc.read_10x_mtx(os.path.join(data_dir, "INF1/filtered_feature_bc_matrix"), gex_only = False)
INF2 = sc.read_10x_mtx(os.path.join(data_dir, "INF2/filtered_feature_bc_matrix"), gex_only = False)

###############################################################################

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


## Slice adata.X to contain only protein data
naive = naive[:, naive.var["feature_types"] == "Gene Expression"].copy()
INF1 = INF1[:, INF1.var["feature_types"] == "Gene Expression"].copy()
INF2 = INF2[:, INF2.var["feature_types"] == "Gene Expression"].copy()

###############################################################################

# BASIC FILTERING

naive.var_names_make_unique()
INF1.var_names_make_unique()
INF2.var_names_make_unique()

sc.pp.filter_cells(naive, min_genes=200)
sc.pp.filter_cells(INF1, min_genes=200)
sc.pp.filter_cells(INF2, min_genes=200)

###############################################################################

# ANNOTATE MITOCHONDRIAL AND RIBOSOMAL GENES
naive.var['mt'] = naive.var_names.str.startswith('mt-')
naive.var['ribosomal'] = naive.var_names.str.match('^(Rpl|Rps)\\d+')
sc.pp.calculate_qc_metrics(naive, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

INF1.var['mt'] = INF1.var_names.str.startswith('mt-')
INF1.var['ribosomal'] = INF1.var_names.str.match('^(Rpl|Rps)\\d+')
sc.pp.calculate_qc_metrics(INF1, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

INF2.var['mt'] = INF2.var_names.str.startswith('mt-')
INF2.var['ribosomal'] = INF2.var_names.str.match('^(Rpl|Rps)\\d+')
sc.pp.calculate_qc_metrics(INF2, qc_vars=['mt', 'ribosomal'], percent_top=None, log1p=False, inplace=True)

###############################################################################

# VISUALIZE QC DATA

adata_list = [naive, INF1, INF2]

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])
    
###############################################################################

# FILTER CELLS BASED ON QC THRESHOLDS

## Total counts

### Floor threshold
for i, adata in enumerate(adata_list):
    print(f"Dataset {i} shape before total_counts floor threshold: {adata.shape}")
    total_counts_floor = adata.obs['total_counts'].quantile(0.05)
    adata_list[i] = adata[adata.obs['total_counts'] > total_counts_floor, :]
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
    ngenes_floor = adata.obs['n_genes_by_counts'].quantile(0.025)
    adata_list[i] = adata[adata.obs['n_genes_by_counts'] > ngenes_floor, :]
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
naive = naive[naive.obs.pct_counts_mt < 5, :].copy()
INF1 = INF1[INF1.obs.pct_counts_mt < 5, :].copy()
INF2 = INF2[INF2.obs.pct_counts_mt < 5, :].copy()

###############################################################################

# RECALCULATE QC METRICS

sc.pp.calculate_qc_metrics(naive, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(INF1, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pp.calculate_qc_metrics(INF2, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

## Print QC summary data
qc_metrics = ['total_counts', 'pct_counts_mt', 'pct_counts_ribosomal', 'n_genes_by_counts']

for i, adata in enumerate(adata_list, start=1):
    print(f"Dataset {i}:")
    for metric in qc_metrics:
        mean_value = adata.obs[metric].mean()
        print(f"Mean {metric}: {mean_value}")
    print("-" * 30) 
    
###############################################################################

# VISUALIZE QC DATA

adata_list = [naive, INF1, INF2]

for adata in adata_list:
    sc.pl.violin(adata, ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal'])
    
###############################################################################

# DOUBLET DETECTION

adata_list = [naive, INF1, INF2]
adata_list

def run_scrublet(adata_list):
    # Initialize an empty list to store DataFrame rows
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

# Run Scrublet analysis
scrublet_result = run_scrublet(adata_list)

###############################################################################

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
naive, INF1, INF2 = adata_list

for i, adata in enumerate(adata_list):
    print(f"Dataset {i + 1} shape after doublet filter: {adata.shape}")

###############################################################################

# CONCATENATION

## Perform concatenation
adata = naive.concatenate(INF1, INF2, batch_key='sample', batch_categories=['naive', 'INF1', 'INF2'] )
adata.obs

## Annotate experimental condition
adata.obs['group'] = adata.obs['sample'].apply(lambda x: 'naive' if x == 'naive' else 'infected')
print(adata.obs[['sample', 'group']].value_counts())


###############################################################################

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

## Select highly variable genes
sc.pp.highly_variable_genes(adata, n_top_genes=3000, subset=True, layer='counts', flavor="seurat_v3")

## Scale data
sc.pp.scale(adata, max_value=10)
sc.pp.scale(bdata, max_value=10)

adata.layers['scaled'] = adata.X.copy()
bdata.layers['scaled'] = bdata.X.copy()

###############################################################################

# Check Nos2

gene_name = "Nos2"
if gene_name in adata.var_names:
    print(f"{gene_name} is present in the dataset.")
else:
    print(f"{gene_name} is not found in the dataset.")


###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/pre-processed"

adata.X = csr_matrix(adata.X)
bdata.X = csr_matrix(bdata.X)

#sc.write(os.path.join(save_dir,'pp-concat-hvg.h5ad'), adata)
#sc.write(os.path.join(save_dir,'pp-concat-full.h5ad'), bdata)

adata.write_h5ad(os.path.join(save_dir, 'pp-concat-hvg.h5ad'))
bdata.write_h5ad(os.path.join(save_dir, 'pp-concat-full.h5ad'))

###############################################################################


