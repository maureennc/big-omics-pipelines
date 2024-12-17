#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 13 08:45:34 2024

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

print(sns.__version__)
print(pd.__version__)
print(np.__version__)
print(sc.__version__)
print(scvi.__version__)

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

## Scanpy
sc.set_figure_params(scanpy = True, dpi = 100, dpi_save = 200, fontsize = 14, figsize = None)

################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/3/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '3-trained-annotated.h5ad'))

################################################################################################################################

# IMPORT MODEL

scvi_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/3/scvi'

model_A_dir = os.path.join(scvi_dir, 'model_A')
print(model_A_dir)
model_A = scvi.model.SCVI.load(model_A_dir, adata=adata)
model_A


################################################################################################################################

# FIND MARKER GENES FOR EACH CLUSTER

sc.tl.rank_genes_groups(adata, groupby='cluster', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 30)
markers = sc.get.rank_genes_groups_df(adata, None)

## Filter to top 50 significant hits
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] 
markers = markers.groupby('group').apply(lambda x: x.nlargest(50, 'scores')).reset_index(drop=True)

save_dir = '/Users/maureen/Documents/projects/harris-lab/maureen/merfish-panel-2/design/spreadsheets/sc-seq'
markers.to_csv(os.path.join(save_dir, 'sig-cluster-markers.csv'), index = False)

################################################################################################################################

# DIFFERENTIAL EXPRESSION - MICROGLIA vs. MACROPHAGES (DURING INFECTION)

scvi_de1 = model_A.differential_expression(idx1=(adata.obs['cell_type'] == 'Macrophage') & (adata.obs['group'] == 'Infected'),
                                          idx2=(adata.obs['cell_type'] == 'Microglia') & (adata.obs['group'] == 'Infected'))

## Filter
scvi_de1 = scvi_de1[(scvi_de1['is_de_fdr_0.05']) & 
                  (abs(scvi_de1.lfc_mean) > 0.5) &
                  ((scvi_de1['raw_normalized_mean1'] + scvi_de1['raw_normalized_mean2']) >= 1)]

scvi_de1.to_csv(os.path.join(save_dir, 'mg-mac-inf-only-de.csv'), index = True)


## Extract list of DE genes and plot
genes = scvi_de1.index.tolist()

os.chdir('/Users/maureen/Documents/projects/harris-lab/maureen/merfish-panel-2/design/figures')

#for gene in genes:
    #sc.pl.umap(adata, color=[gene], save=f'_{gene}-mg-mac-inf.png', show=False)
    
################################################################################################################################

# DIFFERENTIAL EXPRESSION - MICROGLIA vs. MACROPHAGES (NAIVE)

scvi_de2 = model_A.differential_expression(idx1=(adata.obs['cell_type'] == 'Macrophage') & (adata.obs['group'] == 'Control'),
                                          idx2=(adata.obs['cell_type'] == 'Microglia') & (adata.obs['group'] == 'Control'))

## Filter
scvi_de2 = scvi_de2[(scvi_de2['is_de_fdr_0.05']) & 
                  (abs(scvi_de2.lfc_mean) > 0.5) &
                  ((scvi_de2['raw_normalized_mean1'] + scvi_de2['raw_normalized_mean2']) >= 1)]

scvi_de2.to_csv(os.path.join(save_dir, 'mg-mac-naive-only-de.csv'), index = True)

## Extract list of DE genes and plot
genes = scvi_de2.index.tolist()

os.chdir('/Users/maureen/Documents/projects/harris-lab/maureen/merfish-panel-2/design/figures')

#for gene in genes:
    #sc.pl.umap(adata, color=[gene], save=f'_{gene}-mg-mac-naive.png', show=False)

################################################################################################################################

# DIFFERENTIAL EXPRESSION - NAIVE vs. INFECTED MICROGLIA

scvi_de3 = model_A.differential_expression(idx1=(adata.obs['cell_type'] == 'Microglia') & (adata.obs['group'] == 'Infected'),
                                          idx2=(adata.obs['cell_type'] == 'Microglia') & (adata.obs['group'] == 'Control'))

## Filter
scvi_de3 = scvi_de3[(scvi_de3['is_de_fdr_0.05']) & 
                  (abs(scvi_de3.lfc_mean) > 0.5) &
                  ((scvi_de3['raw_normalized_mean1'] + scvi_de3['raw_normalized_mean2']) >= 1)]

scvi_de3.to_csv(os.path.join(save_dir, 'mg-naive-vs-inf-de.csv'), index = True)

## Extract list of DE genes and plot
genes = scvi_de3.index.tolist()

os.chdir('/Users/maureen/Documents/projects/harris-lab/maureen/merfish-panel-2/design/figures')

#for gene in genes:
    #sc.pl.umap(adata, color=[gene], save=f'_{gene}-mg-naive-inf.png', show=False)

################################################################################################################################

# DIFFERENTIAL EXPRESSION - NAIVE vs. INFECTED MACROPHAGE

scvi_de4 = model_A.differential_expression(idx1=(adata.obs['cell_type'] == 'Macrophage') & (adata.obs['group'] == 'Infected'),
                                          idx2=(adata.obs['cell_type'] == 'Macrophage') & (adata.obs['group'] == 'Control'))

## Filter
scvi_de4 = scvi_de4[(scvi_de4['is_de_fdr_0.05']) & 
                  (abs(scvi_de4.lfc_mean) > 0.5) &
                  ((scvi_de4['raw_normalized_mean1'] + scvi_de4['raw_normalized_mean2']) >= 1)]

scvi_de4.to_csv(os.path.join(save_dir, 'mac-naive-vs-inf-de.csv'), index = True)

## Extract list of DE genes and plot
genes = scvi_de4.index.tolist()

os.chdir('/Users/maureen/Documents/projects/harris-lab/maureen/merfish-panel-2/design/figures')

#for gene in genes:
    #sc.pl.umap(adata, color=[gene], save=f'_{gene}-mac-naive-inf.png', show=False)

################################################################################################################################

# COMPILE RESULTS INTO SINGLE DF

scvi_de_combined = pd.concat([scvi_de1, scvi_de2, scvi_de3, scvi_de4])

scvi_de_combined['comparison'] = ['mg-mac-inf'] * len(scvi_de1) + \
                                 ['mg-mac-naive'] * len(scvi_de2) + \
                                 ['mg-naive-inf'] * len(scvi_de3) + \
                                 ['mac-naive-inf'] * len(scvi_de4)

gene_counts = scvi_de_combined.groupby(scvi_de_combined.index).size().reset_index(name='comparison_count')

scvi_de_final = scvi_de_combined.merge(gene_counts, left_index=True, right_on='index')
scvi_de_final = scvi_de_final.sort_values(by='bayes_factor', ascending=False).drop_duplicates(subset='index', keep='first')

scvi_de_final.set_index('index', inplace=True)

scvi_de_final.to_csv(os.path.join(save_dir, 'compiled_de_results.csv'))

################################################################################################################################

# VISUALIZATION

genes = ['Siglech', 'P2ry12', 'Slco2b1', 'Rasgrp3', 'Olfml3', 'Plac8', 'Cx3cr1', 'Zfhx3', 'Agmo', 'Camk1', 'Arhgap22']

sc.pl.umap(adata, color = ['Ifngr1',
'Ifngr2',
'4930430E12Rik',
'4933406I18Rik',
'A630001O12Rik',
'AA467197'])

################################################################################################################################