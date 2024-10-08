#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 14 11:37:00 2024

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

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad/eda"

adata = sc.read_h5ad(os.path.join(data_dir, 'model_H.h5ad'))

###############################################################################

# INITIAL CLUSTERING

sc.pp.neighbors(adata, use_rep = 'X_scVI_H', random_state = 0)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added='leiden_scVI', resolution=1)

sc.pl.umap(adata, color = ['group'])
sc.pl.umap(adata, color = ['leiden_scVI'], legend_loc = 'on data')


###############################################################################

# CHECK CLUSTER MEAN DOUBLET SCORES

## Doublet score by cluster
def print_average_doublet_scores(adata, cluster_col='leiden_scVI', doublet_score_col='doublet_scores'):
    average_doublet_scores = adata.obs.groupby(cluster_col)[doublet_score_col].mean()
    for cluster, score in average_doublet_scores.items():
        print(f"Average doublet score for cluster {cluster}: {score}")

print_average_doublet_scores(adata)

###############################################################################

# GET CLUSTER GENE MARKERS

## Plot
sc.tl.rank_genes_groups(adata, groupby='leiden_scVI', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 30, fontsize = 6)

## Spreadsheet
#markers = sc.get.rank_genes_groups_df(adata, None)
#markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)] # Keep all results
#markers


###############################################################################

# VISUALIZE GENE EXPRESSION

genes = ['Rbfox3', 'Map2', 'Syn1',
         'Slc17a6', 'Slc17a7', 'Nrgn', 'Gad1', 'Gad2',
         'Sox10', 'Pdgfra', 'Cspg4',
         'Col1a1', 'Col1a2', 'Lama1', 
         'Aqp4', 'Gfap', 'Aldh1l1',
         'Hexb', 'Cx3cr1', 'P2ry12',
         'Mbp', 'Plp1', 'Mog',
         'Ttr', 'Gas6', 'Foxj1', 'Aqp1']


sc.pl.dotplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)
sc.pl.matrixplot(adata, genes, groupby = 'leiden_scVI', standard_scale = 'var', use_raw = True, dendrogram = True)

###############################################################################

# ANNOTATIONS

cell_class = { 
"0": "Astrocyte",
"1": "Neuron",
"2": "Neuron",
"3": "Oligodendrocyte",
"4": "Neuron",
"5": "Astrocyte",
"6": "Choroid-plexus epithelial", # Ttr+
"7": "Neuron",
"8": "Neuron",
"9": "Neuron",
"10": "OPC",
"11": "Neuron", 
"12": "Neuron", 
"13": "Neuron", 
"14": "Neuron", 
"15": "Neuron",
"16": "Neuron", 
"17": "Neuron",
"18": "Neuron", 
"19": "Fibroblast",
"20": "Microglia", 
"21": "Neuron", 
"22": "Neuron", 
"23": "Neuron", 
"24": "Neuron",
"25": "Neuron",
"26": "Astrocyte", 
"27": "Neuron", 
"28": "Artifact", # Artifact; high-ribosomal counts in Group F only
"29": "Neuron", 
"30": "Artifact", # Potential artifact cluster; highest doublet score
"31": "Neuron", 
"32": "Unknown",
"33": "Neuron",
"34": "Neuron",
"35": "Artifact", # Potential microglia-oligo doublets
"36": "Artifact" # Ttr-astrocyte doublets
}

adata.obs['cell_class'] = adata.obs.leiden_scVI.map(cell_class)
sc.pl.umap(adata, color = 'cell_class')



cell_type= { 
"0": "Astrocyte",
"1": "Excitatory neuron",
"2": "Excitatory neuron",
"3": "Oligodendrocyte",
"4": "Excitatory neuron",
"5": "Astrocyte",
"6": "Choroid-plexus epithelial", # Ttr+
"7": "Inhibitory neuron",
"8": "Inhibitory neuron",
"9": "Excitatory neuron",
"10": "OPC",
"11": "Inhibitory neuron", 
"12": "Mixed neuron", 
"13": "Excitatory neuron", 
"14": "Excitatory neuron", 
"15": "Excitatory neuron",
"16": "Excitatory neuron", 
"17": "Excitatory neuron",
"18": "Inhibitory neuron", 
"19": "Fibroblast",
"20": "Microglia", 
"21": "Excitatory neuron", 
"22": "Inhibitory neuron", 
"23": "Inhibitory neuron", 
"24": "Excitatory neuron",
"25": "Excitatory neuron",
"26": "Astrocyte", 
"27": "Excitatory neuron", 
"28": "Artifact", # Artifact; high-ribosomal counts in Group F only
"29": "Inhibitory neuron", 
"30": "Artifact", # Potential artifact cluster; highest doublet score
"31": "Excitatory neuron", 
"32": "Unknown",
"33": "Excitatory neuron",
"34": "Excitatory neuron",
"35": "Artifact", # Potential microglia-oligo doublets
"36": "Artifact" # Ttr-astrocyte doublets
}

adata.obs['cell_type'] = adata.obs.leiden_scVI.map(cell_type)
sc.pl.umap(adata, color = 'cell_type')



cluster= { 
"0": "Astrocyte 1",
"1": "Excitatory neuron 1",
"2": "Excitatory neuron 2",
"3": "Oligodendrocyte",
"4": "Excitatory neuron 3",
"5": "Astrocyte 2",
"6": "Choroid-plexus epithelial", # Ttr+
"7": "Inhibitory neuron 1",
"8": "Inhibitory neuron 2",
"9": "Excitatory neuron 4",
"10": "OPC",
"11": "Inhibitory neuron 3", 
"12": "Mixed Ex/In", 
"13": "Excitatory neuron 5", 
"14": "Excitatory neuron 6", 
"15": "Excitatory neuron 7",
"16": "Excitatory neuron 8", 
"17": "Excitatory neuron 9",
"18": "Inhibitory neuron 5", 
"19": "Fibroblast",
"20": "Microglia", 
"21": "Excitatory neuron 10", 
"22": "Inhibitory neuron 6", 
"23": "Inhibitory neuron 7", 
"24": "Excitatory neuron 11",
"25": "Excitatory neuron 12",
"26": "Astrocyte 3", 
"27": "Excitatory neuron 13", 
"28": "Artifact_1_HighRibo", # Artifact; high-ribosomal counts in Group F only
"29": "Inhibitory neuron 8", 
"30": "Artifact_2_HighDoublet", # Potential artifact cluster; highest doublet score
"31": "Excitatory neuron 14", 
"32": "Unknown",
"33": "Excitatory neuron 15",
"34": "Excitatory neuron 16",
"35": "Artifact_3_MixedType", # Potential microglia-oligo doublets
"36": "Artifact_4_MixedType" # Ttr-astrocyte doublets
}
adata.obs['cluster'] = adata.obs.leiden_scVI.map(cluster)
sc.pl.umap(adata, color = 'cluster')


sc.pl.dotplot(adata, genes, groupby = 'cluster', standard_scale = 'var', use_raw = True, dendrogram = True)

###############################################################################

# BASIC CELL COMPOSITION SUMMARY

cell_type_group_numbers = pd.DataFrame(adata.obs.groupby('group')['cell_type'].value_counts().unstack())
cluster_group_numbers = pd.DataFrame(adata.obs.groupby('group')['cluster'].value_counts().unstack())
leiden_group_numbers = pd.DataFrame(adata.obs.groupby('group')['leiden_scVI'].value_counts().unstack())

###############################################################################

# EXPAND GROUP ANNOTATIONS
 
condition = {
    'A': "Sham + AAV-GFP",
    'B': "Sham + AAV-VEGFC",
    'C': "TBI + AAV-GFP ipsilateral",
    'D': "TBI + AAV-GFP contralateral",
    'E': "TBI + AAV-VEGFC ipsilateral",
    'F': "TBI + AAV-VEGFC contralateral"
}

adata.obs['condition'] = adata.obs.group.map(condition)

###############################################################################

# EXPORT

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata.X = csr_matrix(adata.X)

adata.write_h5ad(os.path.join(save_dir, '3-tbi-annotated-hvg.h5ad'))

###############################################################################

# FIGURES

os.chdir("/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures")
         
sc.pl.umap(adata, color = 'cell_class', save = "_cell_class.pdf")
sc.pl.umap(adata, color = 'cell_type', save = "_cell_type.pdf")

###############################################################################