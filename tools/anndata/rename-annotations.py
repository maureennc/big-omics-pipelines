#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 13:19:54 2024

@author: maureen
"""

import os
import scanpy as sc
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '4-tbi-annotated-full.h5ad'))

###############################################################################

# PREPARE ANNDATA

## Mask artifact clusters
adata = adata[adata.obs['cell_type'] != 'Unknown'].copy()

## Rename Ttr+
adata.obs['cell_type'].replace({'Ttr+': 'Choroid-plexus epithelial'}, inplace = True)
adata.obs['cluster'].replace({'Ttr+': 'Choroid-plexus epithelial'}, inplace = True)


## Rename Mixed population
adata.obs['cluster'].replace({'Inhibitory neuron 4': 'Mixed Ex/In'}, inplace = True)

## Update cell_type assignment
adata.obs['cell_type'] = adata.obs['cell_type'].cat.add_categories('Mixed neuron')
adata.obs.loc[adata.obs['cluster'] == 'Mixed Ex/In', 'cell_type'] = 'Mixed neuron'

sc.pl.umap(adata, color = 'cell_type', legend_loc = 'on data')
sc.pl.umap(adata, color = 'cluster',  legend_loc = 'on data', legend_fontsize = 5)
sc.pl.umap(adata, color = 'leiden_scVI', legend_loc = 'on data')

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/h5ad"

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, '6-tbi-annotated-full.h5ad'))

###############################################################################