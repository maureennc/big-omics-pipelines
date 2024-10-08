#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 12:39:54 2024

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

###############################################################################

# SETTINGS

## Random seed
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)
scvi.settings.seed = 0

## Matplotlib
%matplotlib qt5
#plt.rcParams['font.family'] = 'Arial'

## Scanpy
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)

###############################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad/eda"
scvi_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/rivanna/scvi"

adata = sc.read_h5ad(os.path.join(data_dir, 'model_H.h5ad'))

###############################################################################

# EXPLORE GENES OF INTEREST

#umap_save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/figures/model_H-annotation"
os.chdir("/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/figures/model_H-annotation")


genes = ['Rbfox3', 'Map2', 'Syn1', 'Slc17a6', 'Slc17a7', 'Gad1', 'Gad2', # Neurons
              'Cx3cr1', 'Hexb', 'P2ry12', 'Mrc1', # Microglia / PVM
              'Gfap', 'Aldh1l1', 'Aqp4', 'Atp1a2', # Astrocytes
              'Pecam1', 'Tie1', 'Cldn5', # Endothelial cells
              'Pdgfrb', 'Rgs5', # Pericytes
              'Pdgfra', 'Cspg4', 'Sox10', # OPC/Oligos
              'Mbp', 'Plp1', 'Mog', 'Itgam', 'Itgax', # Oligos / Microglia
              'Htr2c', 'Enpp2', 'Stk39', 'Ttr', 'Cdh1', 'Aqp1', 'Foxj1', 'Krt18', 'Gas6', # CPE
              'Vim', 'Col1a2', 'Col1a1', 'Col3a1', 'Acta2', 'Myh11', 'Thy1', 'Col4a1', 'Lama1', 'Lama2', 'Fn1'] # Fibroblasts

for gene in genes:
    sc.pl.umap(adata, color=gene, save= f'_{gene}.pdf')
    
    
params = ['total_counts', 'n_genes_by_counts', 'pct_counts_mt', 'pct_counts_ribosomal', 'doublet_scores']

for param in params:
    sc.pl.umap(adata, color=param, save= f'_{param}.pdf')

sc.pl.umap(adata, color = 'leiden_scVI', legend_loc = 'on data', legend_fontsize = 10, save = '_leiden_legend_on-data.pdf')

###############################################################################

# GROUP UMAPs

sc.pl.umap(adata, color = 'condition', save = "_condition.pdf")

sc.pl.umap(adata, color = 'condition', groups = ['Sham + AAV-GFP'], save = "_treatment_condition_A.pdf")
sc.pl.umap(adata, color = 'condition', groups = ['Sham + AAV-VEGFC'], save = "_treatment_condition_B.pdf")
sc.pl.umap(adata, color = 'condition', groups = ['TBI + AAV-GFP ipsilateral'], save = "_treatment_condition_C.pdf")
sc.pl.umap(adata, color = 'condition', groups = ['TBI + AAV-GFP contralateral'], save = "_treatment_condition_D.pdf")
sc.pl.umap(adata, color = 'condition', groups = ['TBI + AAV-VEGFC ipsilateral'], save = "_treatment_condition_E.pdf")
sc.pl.umap(adata, color = 'condition', groups = ['TBI + AAV-VEGFC contralateral'], save = "_treatment_condition_F.pdf")


###############################################################################

# PLOT UPDATED DOUBLET SUMMARY

sc.pl.violin(adata, 'doublet_scores', groupby = 'group', save = "_doublet_scores.pdf")


data_for_plot = []

for group_label in adata.obs['group'].unique():
    group_size = adata.obs['group'].value_counts()[group_label]
    data_for_plot.append([group_label, group_size])


df = pd.DataFrame(data_for_plot, columns=['Dataset', 'Cells'])

## Bar graph
fig, ax = plt.subplots()
df.plot(x='Dataset', y='Cells', kind='bar', ax=ax, color='blue')
ax.set_title('Treatment group counts')
ax.set_xlabel('Group')
ax.set_ylabel('Number of cells')
plt.xticks(rotation=0)
plt.grid(False)
plt.legend().set_visible(False)  # Hide the legend
plt.tight_layout()
plt.savefig("number_cells_per_treatment_group.pdf", format='pdf')
plt.show()

###############################################################################