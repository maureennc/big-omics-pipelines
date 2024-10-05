#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 09:15:21 2024

@author: maureen
"""

import os
import scanpy as sc
import scvi
import pandas as pd

###############################################################################

# DEFINE PATHS

input_dir = '/scratch/mnc3ra/tbi_snseq/solo/input'
output_dir = '/scratch/mnc3ra/tbi_snseq/solo/output'

input_file = '1-raw-concat.h5ad'
output_file = 'solo-results-full.csv'

###############################################################################

# IMPORT DATA

adata = sc.read_h5ad(os.path.join(input_dir, input_file))

###############################################################################

# PREPARE DATA

sc.pp.filter_genes(adata, min_cells = 10)
sc.pp.highly_variable_genes(adata, n_top_genes = 3000, subset = True, flavor = 'seurat_v3')

###############################################################################

# TRAIN VAE

## VAE
scvi.model.SCVI.setup_anndata(adata)
vae = scvi.model.SCVI(adata)
vae.train()

###############################################################################

## TRAIN AND RUN SOLO

all_results = pd.DataFrame()

for doublet_ratio in [1, 2, 3]:
    
    ## Train model
    solo = scvi.external.SOLO.from_scvi_model(vae, adata=adata, doublet_ratio=doublet_ratio)
    solo.train()

    ## Extract predictions
    results = solo.predict()
    results['prediction_soft'] = solo.predict(soft=False)
    results['doublet_ratio'] = doublet_ratio
    results['group'] = results.index.str.split('-').str[-1]

    all_results = pd.concat([all_results, results], ignore_index=True)


## Export results
all_results.to_csv(os.path.join(output_dir, output_file), index=True)

###############################################################################