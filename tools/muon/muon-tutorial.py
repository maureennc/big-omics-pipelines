#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar  9 11:31:54 2024

@author: maureen
"""

# Change directory to the root folder of the repository
import os
os.chdir("../")

# Download data

# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3

# This is the directory where those files are downloaded to
data_dir = "data/pbmc5k_protein"

for file in os.listdir(data_dir):
    print(file)

###############################################################################

# LOAD LIBRARIES AND DATA

import numpy as np
import pandas as pd
import scanpy as sc
import muon as mu
from muon import prot as pt

mdata = mu.read_10x_mtx(os.path.join(data_dir, "filtered_feature_bc_matrix"))
mdata_raw = mu.read_10x_mtx(os.path.join(data_dir, "raw_feature_bc_matrix"))

###############################################################################

# PROTEIN

prot = mdata.mod['prot']
prot

###############################################################################

# DSB NORMALIZATION

pt.pp.dsb(mdata, raw=mdata_raw, empty_droplets=droplets)
mdata_raw['rna'].obs["log10umi"] = np.array(np.log10(mdata_raw['rna'].X.sum(axis=1) + 1)).reshape(-1)
mu.pl.histogram(mdata_raw['rna'], ['log10umi'], bins=50)
mu.pl.histogram(mdata_raw['rna'][mdata_raw['rna'].obs.log10umi >= 1], ['log10umi'], bins=50)

isotypes = mdata_raw['prot'].var_names[29:32].values
isotypes

prot.layers['counts'] = prot.X

pt.pp.dsb(mdata, mdata_raw, empty_counts_range=(1.5, 2.8), isotype_controls=isotypes, random_state=1)

sc.pl.scatter(mdata['prot'], x="CD3_TotalSeqB", y="CD19_TotalSeqB", layers='counts')
sc.pl.scatter(mdata['prot'], x="CD3_TotalSeqB", y="CD19_TotalSeqB")

###############################################################################

# DOWNSTREAM ANALYSIS

sc.tl.pca(prot)
sc.pl.pca(prot, color='CD3_TotalSeqB')

sc.pp.neighbors(prot)
sc.tl.umap(prot, random_state=1)

sc.pl.umap(prot, color=['CD3_TotalSeqB', 'CD14_TotalSeqB'])

###############################################################################

# RNA

rna = mdata.mod['rna']
rna

## QC
rna.var['mt'] = rna.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

sc.pl.violin(rna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)


## Filter genes which expression is not detected:
mu.pp.filter_var(rna, 'n_cells_by_counts', lambda x: x >= 3)
# Same as the following but doesn't copy the object:
#   sc.pp.filter_genes(rna, min_cells=3)

mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: (x >= 200) & (x < 5000))
# Same as the following but doesn't copy the object
#   sc.pp.filter_cells(rna, min_genes=200)
#   rna = rna[rna.obs.n_genes_by_counts < 5000, :]


## Filter cells
mu.pp.filter_obs(rna, 'total_counts', lambda x: (x > 1500) & (x < 15000))
mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < 20)

## Visualize QC data
sc.pl.violin(rna, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True)

###############################################################################

# Normalization
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

# Feature selection

sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)
sc.pl.highly_variable_genes(rna)

# Scaling
rna.raw = rna
sc.pp.scale(rna, max_value=10)

###############################################################################


