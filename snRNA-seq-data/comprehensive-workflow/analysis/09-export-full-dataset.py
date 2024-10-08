#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:28:42 2024

@author: maureen
"""

import os
import pandas as pd
import scanpy as sc
import scipy.sparse
from scipy.sparse import csr_matrix

################################################################################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '4-tbi-annotated-full.h5ad'))
adata.X = adata.layers['log1p'].copy()

################################################################################################################################

# PREPARE DATA

## Exclude artifacts
adata = adata[~adata.obs['cell_type'].isin(['Artifact', 'Unknown'])].copy()

## Update annotations
adata.obs['cell_type'] = adata.obs['cell_type'].replace('Ttr+', 'Unknown')
adata.obs['cluster'] = adata.obs['cluster'].replace('Ttr+', 'Unknown')

## Alphabetize annotations
adata.obs['cluster'] = pd.Categorical(adata.obs['cluster'], categories=sorted(adata.obs['cluster'].unique()), ordered=True) # commented 6/25/24
#adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'], categories=sorted(adata.obs['cell_type'].unique()), ordered=True)

## Renumber leiden scores
value_counts = adata.obs['leiden_scVI'].value_counts()

cluster_mapping = {old_label: rank for rank, old_label in enumerate(value_counts.index)}
adata.obs['leiden'] = adata.obs['leiden_scVI'].map(cluster_mapping)
print(adata.obs['leiden'].value_counts())

################################################################################################################################

# EXPORT COUNTS MATRIX LAYERS

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/csv/full/data'

def save_matrix(matrix, obs_names, var_names, file_path):
    if isinstance(matrix, scipy.sparse.spmatrix):
        df = pd.DataFrame.sparse.from_spmatrix(matrix, index=obs_names, columns=var_names)
    else:
        df = pd.DataFrame(matrix, index=obs_names, columns=var_names)
    df.to_csv(file_path, sep=',')


layers_to_export = ['counts', 'log1p', 'normalized']

for layer_name in layers_to_export:
    if layer_name in adata.layers:
        layer_matrix = adata.layers[layer_name]
        layer_file_path = os.path.join(save_dir, f'{layer_name}-matrix.csv')
        save_matrix(layer_matrix, adata.obs_names, adata.var_names, layer_file_path)
        print(f"Exported {layer_name} layer")
    else:
        print(f"'{layer_name}' layer not saved in AnnData")


################################################################################################################################

# EXPORT METADATA

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/csv/full/metadata'

adata.obs.to_csv(os.path.join(save_dir, 'cell-metadata.csv'))
adata.var.to_csv(os.path.join(save_dir, 'feature-metadata.csv'))

################################################################################################################################

# EXPORT DIMENSINONAL REDUCTIONS

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/csv/full/reductions'

pca = pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names)
pca.to_csv(os.path.join(save_dir, 'pca.csv'))

latent_representation = pd.DataFrame(adata.obsm['X_scVI_H'], index=adata.obs_names)
latent_representation.to_csv(os.path.join(save_dir, 'latent_representation.csv'))


umap = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs_names)
umap.to_csv(os.path.join(save_dir, 'umap.csv'))


################################################################################################################################

# NEIGHBORS GRAPH
## Takes a long time; can recalculate in Seurat

#connectivities = pd.DataFrame.sparse.from_spmatrix(adata.obsp['connectivities'], index=adata.obs_names, columns=adata.obs_names)
#connectivities.to_csv(os.path.join(save_dir, 'connectivities.csv'))

#distances = pd.DataFrame.sparse.from_spmatrix(adata.obsp['distances'], index=adata.obs_names, columns=adata.obs_names)
#distances.to_csv(os.path.join(save_dir, 'distances.csv'))

################################################################################################################################

# EXPORT UPDATED H5AD

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/h5ad'

adata.X = csr_matrix(adata.X)
adata.write_h5ad(os.path.join(save_dir, '2024_tbi-snrna-seq-cleaned-full.h5ad'))

################################################################################################################################
