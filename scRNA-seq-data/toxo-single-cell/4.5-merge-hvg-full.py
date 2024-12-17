#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 15:01:16 2024

@author: maureen
"""

import os
import scanpy as sc
from scipy.sparse import csr_matrix

###############################################################################

# IMPORT

data_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/annotated-trained"
adata = sc.read_h5ad(os.path.join(data_dir, "pp-annotated-model_A.h5ad"))

data_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/pre-processed'
bdata = sc.read_h5ad(os.path.join(data_dir, "pp-concat-full.h5ad"))


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

data_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/pre-processed'

bdata.X = csr_matrix(bdata.X)

bdata.write_h5ad(os.path.join(data_dir, 'pp-concat-full-metadata.h5ad'))

###############################################################################