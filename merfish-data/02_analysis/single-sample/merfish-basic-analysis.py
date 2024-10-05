#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 17:47:47 2024

@author: maureen
"""

import os
import scanpy as sc
from scanpy import tl, pl, pp

import squidpy as sq
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 200

################################################################################

# IMPORT

data_dir = os.path.join('Documents', 'spatial', 'data', 'vpt', 'E003-v2')

os.listdir(data_dir)

adata = sq.read.vizgen(data_dir,
                       counts_file = "E003-v2_cell_by_gene.csv",
                       meta_file = 'E003-v2_cell_metadata.csv',
                       transformation_file = 'micron_to_mosaic_pixel_transform.csv')

genes = adata.var
genes

################################################################################

# QC

## Filter blanks
mask = adata.var_names.str.startswith('Blank')
adata = adata[:, ~mask].copy()

## Calculate QC metrics
sc.pp.calculate_qc_metrics(adata, percent_top = None, log1p = False, inplace = True)
bdata = adata.copy()

## Total counts
sc.pl.violin(adata, 'total_counts')
adata = adata[adata.obs['total_counts'] < 300].copy()
adata = adata[adata.obs['total_counts'] > 50].copy()
sc.pl.violin(adata, 'total_counts')

## Number of genes per cell
sc.pl.violin(adata, 'n_genes_by_counts')
adata = adata[adata.obs['n_genes_by_counts'] < 90].copy()
adata = adata[adata.obs['n_genes_by_counts'] > 20].copy()
sc.pl.violin(adata, 'n_genes_by_counts')

## Volume
sc.pl.violin(adata, 'volume')
adata = adata[adata.obs['volume'] < 1250].copy()
adata = adata[adata.obs['volume'] > 250].copy()
sc.pl.violin(adata, 'volume')

################################################################################

# NORMALIZATION AND VARIANCE STABILIZATION

adata.layers['counts'] = adata.X.copy()

## Normalization
target_sum = adata.obs['total_counts'].median()
sc.pp.normalize_total(adata, target_sum = target_sum)
adata.layers['normalized'] = adata.X.copy()

## Log1p
sc.pp.log1p(adata)
adata.layers['log1p'] = adata.X.copy()
adata.raw = adata

## Scale
sc.pp.scale(adata)
adata.layers['scaled'] = adata.X.copy()

adata.X = adata.layers['log1p'].copy()

################################################################################

# CLUSTERING

#adata.X = adata.layers['counts'].copy()
#sc.pp.highly_variable_genes(adata, flavor = 'seurat_v3')

sc.pp.pca(adata)
sc.pl.pca_variance_ratio(adata)
sc.pp.neighbors(adata, n_pcs = 10)
sc.tl.leiden(adata, resolution = 1)
sc.tl.umap(adata)

sc.pl.umap(adata, color = 'leiden')


################################################################################

# ANNOTATION

cell_type = {str(i): "" for i in range(21)}

# annotated dictionary here

#adata.obs['cell_type'] = adata.obs.['leiden'].map(cell_type)

################################################################################

# CLUSTER DIFFERENTIAL EXPRESSION

groups = [str(i) for i in range(11)]

subset = adata[adata.obs['leiden'].isin(groups)].copy()

sc.tl.rank_genes_groups(adata, groupby = 'leiden')
markers = sc.get.rank_genes_groups_df(adata, group = None)
markers

top_markers = markers.groupby('group').head(5)
top_markers = top_markers['names'].tolist()
sc.pl.heatmap(adata, top_markers, groupby = 'leiden')

################################################################################

# SPATIAL OVERVIEW

sq.gr.spatial_neighbors(adata)

sq.pl.spatial_scatter(adata, color = 'leiden', shape = None, groups = ['0', '1', '2'])

sq.gr.spatial_autocorr(adata)
adata

moran = pd.DataFrame(adata.uns['moranI'])

moran
sq.pl.spatial_scatter(adata, shape = None, color = 'Slc17a7', cmap = 'rocket')


sq.gr.nhood_enrichment(adata, cluster_key = 'leiden')
sq.pl.nhood_enrichment(adata, cluster_key = 'leiden', figsize = (3,3))

markers.set_index('group', inplace = True)
markers.loc['17']

################################################################################