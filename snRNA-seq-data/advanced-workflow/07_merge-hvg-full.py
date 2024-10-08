#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 15:15:19 2024

@author: maureen
"""

import os
import scanpy as sc
from scipy.sparse import csr_matrix

###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '3-tbi-annotated-hvg.h5ad'))
bdata = sc.read_h5ad(os.path.join(data_dir, '2-tbi-seq-full.h5ad'))

###############################################################################

# FIND INTERSECTION OF BARCODES

common_barcodes = adata.obs_names.intersection(bdata.obs_names)
bdata = bdata[common_barcodes].copy()

###############################################################################

# TRANSFER METADATA TO GENOME-SCALE BRANCH

## adata.obs
unique_columns = [col for col in adata.obs.columns if col not in bdata.obs.columns]

bdata.obs = bdata.obs.join(adata.obs[unique_columns], how='left')

## uns entries
for key in adata.uns.keys():
    if key not in bdata.uns:
        bdata.uns[key] = adata.uns[key]

## obsm entries
for key in adata.obsm.keys():
    if key not in bdata.obsm:
        bdata.obsm[key] = adata.obsm[key]

## obsp entries
for key in adata.obsp.keys():
    if key not in bdata.obsp:
        bdata.obsp[key] = adata.obsp[key]

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

bdata.X = csr_matrix(bdata.X)

bdata.write_h5ad(os.path.join(save_dir, '4-tbi-annotated-full.h5ad'))

###############################################################################