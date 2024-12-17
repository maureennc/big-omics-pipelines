#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:16:21 2024

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

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/3/h5ad"

## adata used for doublet training
adata_naive = sc.read_h5ad(os.path.join(data_dir, '1-naive-qc-filtered.h5ad'))
adata_INF1 = sc.read_h5ad(os.path.join(data_dir, '1-INF1-qc-filtered.h5ad'))
adata_INF2 = sc.read_h5ad(os.path.join(data_dir, '1-INF2-qc-filtered.h5ad'))


## bdata used for processing and exports
bdata_naive = adata_naive.copy()
bdata_INF1 = adata_INF1.copy()
bdata_INF2 = adata_INF2.copy()

################################################################################################################################

# PREPARE DATA FOR TRAINING

## Dimensionality reduction
sc.pp.highly_variable_genes(adata_naive, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')
sc.pp.highly_variable_genes(adata_INF1, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')
sc.pp.highly_variable_genes(adata_INF2, n_top_genes = 2000, subset = True, flavor = 'seurat_v3')

## Made adata.X CSR
adata_naive.X = csr_matrix(adata_naive.X)
adata_INF1.X = csr_matrix(adata_INF1.X)
adata_INF2.X = csr_matrix(adata_INF2.X)

################################################################################################################################

# DOUBLET DETECTION - Naive

## Naive - Train VAE model
scvi.model.SCVI.setup_anndata(adata_naive)
vae_naive = scvi.model.SCVI(adata_naive)
vae_naive.train()

## Naive - Train SOLO model
solo_naive = scvi.external.SOLO.from_scvi_model(scvi_model = vae_naive,
                                                adata = adata_naive,
                                                doublet_ratio = 2) # default doublet_ratio
solo_naive.train()

## Naive - Extract results
solo_naive_results = solo_naive.predict()
solo_naive_results['prediction'] = solo_naive.predict(soft=False)

print(solo_naive_results.prediction.value_counts())

################################################################################################################################

# DOUBLET DETECTION - INF1

## INF1 - Train VAE model
scvi.model.SCVI.setup_anndata(adata_INF1)
vae_INF1 = scvi.model.SCVI(adata_INF1)
vae_INF1.train()

## INF1 - Train SOLO model
solo_INF1 = scvi.external.SOLO.from_scvi_model(scvi_model = vae_INF1,
                                               adata = adata_INF1,
                                               doublet_ratio = 3) # increased due to macrophage-T cell stickiness
solo_INF1.train()

## INF1 - Extract results
solo_INF1_results = solo_INF1.predict()
solo_INF1_results['prediction'] = solo_INF1.predict(soft=False)

print(solo_INF1_results.prediction.value_counts())

################################################################################################################################

# DOUBLET DETECTION - INF2

## INF2 - Train VAE model
scvi.model.SCVI.setup_anndata(adata_INF2)
vae_INF2 = scvi.model.SCVI(adata_INF2)
vae_INF2.train()

## INF2 - Train SOLO model
solo_INF2 = scvi.external.SOLO.from_scvi_model(scvi_model = vae_INF2,
                                               adata = adata_INF2,
                                               doublet_ratio = 4) # increased due to macrophage-T cell stickiness & high lane loading density
solo_INF2.train()

## INF2 - Extract results
solo_INF2_results = solo_INF2.predict()
solo_INF2_results['prediction'] = solo_INF2.predict(soft=False)

print(solo_INF2_results.prediction.value_counts())


################################################################################################################################

# MAP DOUBLETS TO BDATA

## Add annotations from results table to the corresponding adata
adata_naive.obs[['doublet', 'singlet', 'doublet_predictions']] = solo_naive_results[['doublet', 'singlet', 'prediction']]
adata_INF1.obs[['doublet', 'singlet', 'doublet_predictions']] = solo_INF1_results[['doublet', 'singlet', 'prediction']]
adata_INF2.obs[['doublet', 'singlet', 'doublet_predictions']] = solo_INF2_results[['doublet', 'singlet', 'prediction']]

## Define function to map metrics
def map_metrics(source_adata, target_adata):    
    for column in ['doublet', 'singlet', 'doublet_predictions']:
        target_adata.obs[column] = 'unknown'  # Ensure that target adata has columns initialized

    ## Find intersection of indices to map common barcodes
    common_barcodes = target_adata.obs_names.intersection(source_adata.obs_names)
    target_adata.obs.loc[common_barcodes, ['doublet', 'singlet', 'doublet_predictions']] = source_adata.obs.loc[common_barcodes, ['doublet', 'singlet', 'doublet_predictions']]

## Apply mapping for each dataset
map_metrics(adata_naive, bdata_naive)
map_metrics(adata_INF1, bdata_INF1)
map_metrics(adata_INF2, bdata_INF2)

## Check the results
print(bdata_naive.obs[['doublet', 'singlet', 'doublet_predictions']].head())
print(bdata_INF1.obs[['doublet', 'singlet', 'doublet_predictions']].head())
print(bdata_INF2.obs[['doublet', 'singlet', 'doublet_predictions']].head())

################################################################################################################################


# VISUALIZE DOUBLET PREDICTIONS

## Naive
sc.pp.normalize_total(adata_naive)
sc.pp.log1p(adata_naive)
sc.tl.pca(adata_naive, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_naive, log=True, n_pcs = 50)
sc.pp.neighbors(adata_naive, n_pcs = 30)
sc.tl.umap(adata_naive)
sc.tl.leiden(adata_naive, resolution = 0.5)

sc.pl.umap(adata_naive, color = ['leiden'], legend_loc = 'on data')
sc.pl.umap(adata_naive, color = ['doublet_predictions'])


## INF1
sc.pp.normalize_total(adata_INF1)
sc.pp.log1p(adata_INF1)
sc.tl.pca(adata_INF1, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_INF1, log=True, n_pcs = 50)
sc.pp.neighbors(adata_INF1, n_pcs = 30)
sc.tl.umap(adata_INF1)
sc.tl.leiden(adata_INF1, resolution = 0.5)

sc.pl.umap(adata_INF1, color = ['leiden'], legend_loc = 'on data')
sc.pl.umap(adata_INF1, color = ['doublet_predictions'])


## INF2
sc.pp.normalize_total(adata_INF2)
sc.pp.log1p(adata_INF2)
sc.tl.pca(adata_INF2, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata_INF2, log=True, n_pcs = 50)
sc.pp.neighbors(adata_INF2, n_pcs = 30)
sc.tl.umap(adata_INF2)
sc.tl.leiden(adata_INF2, resolution = 0.5)

sc.pl.umap(adata_INF2, color = ['leiden'], legend_loc = 'on data')
sc.pl.umap(adata_INF2, color = ['doublet_predictions'])

## Inspect INF2 cluster 3 (confirmed likely doublets)
sc.tl.rank_genes_groups(adata_INF2, 'leiden', method='wilcoxon')
top_genes_cluster_3 = adata_INF2.uns['rank_genes_groups']['names']['3'][:10] 
sc.pl.dotplot(adata_INF2, var_names=top_genes_cluster_3, groupby='leiden', dendrogram=True)
sc.pl.rank_genes_groups_dotplot(adata_INF2, n_genes=10)

genes = ['Cd3e', 'Cd4', 'Cd8b1', 'C1qa', 'Trem2', 'Cybb', 'Tnfaip2', 'Cx3cr1']
sc.pl.dotplot(adata_INF2, genes, groupby = 'leiden', dendrogram = True)

sc.pl.umap(adata_INF2, color = ['n_genes_by_counts', 'total_counts'])

################################################################################################################################

# FILTER DOUBLETS

## Before filtering - Print median values of n_genes_by_counts and total_counts
print("Before filtering doublets:")
print("Naive - Median n_genes_by_counts:", bdata_naive.obs['n_genes_by_counts'].median())
print("Naive - Median total_counts:", bdata_naive.obs['total_counts'].median())
print("INF1 - Median n_genes_by_counts:", bdata_INF1.obs['n_genes_by_counts'].median())
print("INF1 - Median total_counts:", bdata_INF1.obs['total_counts'].median())
print("INF2 - Median n_genes_by_counts:", bdata_INF2.obs['n_genes_by_counts'].median())
print("INF2 - Median total_counts:", bdata_INF2.obs['total_counts'].median())

## Filter out doublets
bdata_naive = bdata_naive[bdata_naive.obs['doublet_predictions'] == 'singlet'].copy()
bdata_INF1 = bdata_INF1[bdata_INF1.obs['doublet_predictions'] == 'singlet'].copy()
bdata_INF2 = bdata_INF2[bdata_INF2.obs['doublet_predictions'] == 'singlet'].copy()

## After filtering - Print effect of filtering
print("\nAfter filtering doublets:")
print("Naive - Median n_genes_by_counts:", bdata_naive.obs['n_genes_by_counts'].median())
print("Naive - Median total_counts:", bdata_naive.obs['total_counts'].median())
print("INF1 - Median n_genes_by_counts:", bdata_INF1.obs['n_genes_by_counts'].median())
print("INF1 - Median total_counts:", bdata_INF1.obs['total_counts'].median())
print("INF2 - Median n_genes_by_counts:", bdata_INF2.obs['n_genes_by_counts'].median())
print("INF2 - Median total_counts:", bdata_INF2.obs['total_counts'].median())

################################################################################################################################

# CONCATENATION

## Annotate groups
bdata_naive.obs['group'] = 'Control'
bdata_INF1.obs['group'] = 'Infected'
bdata_INF2.obs['group'] = 'Infected'

## Perform concatenation
bdata = bdata_naive.concatenate(bdata_INF1, bdata_INF2, 
                                 batch_key='sample', 
                                 batch_categories=['Control', 'Infected1', 'Infected2'])

bdata.obs

################################################################################################################################

# FINISH PRE-PROCESSING

## Filter genes appearing in few cells
sc.pp.filter_genes(bdata, min_cells=3)

## Finish pre-processing (Normalization, log1p-transformation)
bdata.layers['counts'] = bdata.X.copy()
sc.pp.normalize_total(bdata)
bdata.layers['normalized'] = bdata.X.copy()
sc.pp.log1p(bdata)
bdata.layers['log1p'] = bdata.X.copy()
bdata.raw = bdata.copy()

################################################################################################################################

# SELECT HVGs

bdata_subset = bdata.copy()
sc.pp.highly_variable_genes(bdata_subset, n_top_genes = 5000, subset = True, layer = 'counts', batch_key = 'sample', flavor = "seurat_v3") # raw counts for seurat-v3

sc.pp.scale(bdata, max_value=10)
bdata.layers['scaled'] = bdata.X.copy() # save

## Revert to log1p
bdata.X = bdata.layers['log1p'].copy()

################################################################################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/3/h5ad"

bdata_subset.X = csr_matrix(bdata_subset.X)
bdata.X = csr_matrix(bdata.X)

bdata_subset.obs['doublet'] = bdata_subset.obs['doublet'].astype(float)
bdata_subset.obs['singlet'] = bdata_subset.obs['singlet'].astype(float)

bdata.obs['doublet'] = bdata.obs['doublet'].astype(float)
bdata.obs['singlet'] = bdata.obs['singlet'].astype(float)

bdata_subset.write_h5ad(os.path.join(save_dir, '2-concat-pp-hvg.h5ad'))
bdata.write_h5ad(os.path.join(save_dir, '2-concat-pp-full.h5ad'))

################################################################################################################################
