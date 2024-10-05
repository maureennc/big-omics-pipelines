#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 16:36:20 2023

@author: maureen
"""
import os
os.chdir("/Users/maureen/Documents/projects/lukens-lab/nick")

import scanpy as sc
import squidpy as sq
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from pathlib import Path

# Make sure you have the right file structure
#!mkdir demo_data
#!mkdir tutorial_data/vizgen_data ## This is where the cell_by_gene and cell_metadata files need to be
#!mkdir tutorial_data/vizgen_data/images ## This is where the transformation file needs to be

# Import data and construct adata object
vizgen_dir = Path().resolve() / "demo_data" / "vizgen_data"

adata = sq.read.vizgen(path = vizgen_dir, 
                       counts_file = "cell_by_gene.csv",
                       meta_file = "cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("mt-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"],percent_top= None, inplace=True)

# Calculate QC Metrics
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.distplot(
    adata.obs["total_counts"],
    kde=False,
    ax=axs[0],
)
sns.distplot(
    adata.obs["total_counts"][adata.obs["total_counts"] < 10000],
    kde=False,
    bins=40,
    ax=axs[1],
)
sns.distplot(
    adata.obs["n_genes_by_counts"],
    kde=False,
    bins=60,
    ax=axs[2],
)
sns.distplot(
    adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000],
    kde=False,
    bins=60,
    ax=axs[3],
)

# Filter out cells and genes with minimum counts
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.filter_genes(adata, min_cells=10)


# Re-visualize filtered results

# Annotate highly variable genes (requires scikit-misc)
adata.layers["counts"] = adata.X.copy()
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=4000)
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)


# Visualize annotation on UMAP and spatial coodinates
sc.pl.umap(
    adata,
    color=[
        "total_counts",
        "n_genes_by_counts",
        "leiden",
    ],
    wspace=0.4,
)

sq.pl.spatial_scatter(adata, shape=None,color=[
        "leiden",
    ],
    wspace=0.4,)

sq.pl.spatial_scatter(adata, shape = None, color = ['Cell type'])


# Identify clusters
sc.pl.umap(adata, color = ['leiden'])
sc.tl.rank_genes_groups(adata, 'leiden')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)
markers = sc.get.rank_genes_groups_df(adata, None)
markers = markers[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > .5)]

# Obtain data frame of panel genes
panel = adata.var_names
panel = pd.Index.to_frame(panel)
panel.to_csv('MERFISH_demo_genes.csv', index = False)

# Astrocytes
sc.pl.umap(adata, color = ['leiden', 'Aqp4'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Aldh1l1'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Cxcl14'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Creb3l1'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Mlc1'], legend_loc = "on data")


# Excitatory Neurons
sc.pl.umap(adata, color = ['leiden', 'Slc17a6'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Slc17a7'], legend_loc = "on data")

# Inhibitory Neurons
sc.pl.umap(adata, color = ['leiden', 'Gad1'], legend_loc = "on data")

# Cholinergic Neurons
sc.pl.umap(adata, color = ['leiden', 'Chat'], legend_loc = "on data")

# OPC
sc.pl.umap(adata, color = ['leiden', 'Pdgfra'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Traf4'], legend_loc = "on data")


# Oligodendrocytes
sc.pl.umap(adata, color = ['leiden', 'Ermn'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Gjc3'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Lpar1'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Ndrg1'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Opalin'], legend_loc = "on data")

# Endothelial
sc.pl.umap(adata, color = ['leiden', 'Igf1r'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Fn1'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Rgs5'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Slco1a4'], legend_loc = "on data")

# Pericytes?
sc.pl.umap(adata, color = ['leiden', 'Ace2'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Adora2a'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Rgs5'], legend_loc = "on data")

# Microglia/Macrophage
sc.pl.umap(adata, color = ['leiden', 'Rgs2'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Man1a'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Selplg'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Slc15a3'], legend_loc = "on data")


# Neuron
sc.pl.umap(adata, color = ['leiden', 'Bdnf'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Cbln1'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Cbln2'], legend_loc = "on data")
sc.pl.umap(adata, color = ['leiden', 'Nos1'], legend_loc = "on data")

sc.pl.umap(adata, color = ['leiden', 'Pitx3'], legend_loc = "on data")




## Create a dictionary to map cell labels
Cell_type = {
"0":"Oligodendrocyte",
"1":"Astrocyte",
"2":"Glutamatergic neuron",
"3":"Endothelial cell",
"4":"Glutamatergic neuron",
"5":"Microglia/Macrophage",
"6":"GABAergic neuron",
"7":"Glutamatergic neuron",
"8":"Glutamatergic neuron",
"9":"Glutamatergic neuron",
"10":"OPC",
"11":"Glutamatergic neuron",
"12":"Other neuron",
"13":"GABAergic neuron",
"14":"Astrocyte",
"15":"GABAergic neuron",
"16":"Other neuron",
"17":"Glutamatergic neuron",
"18":"GABAergic neuron",
"19": "Pericyte",
"20": "GABAergic neuron",
"21": "Glutamatergic neuron",
"22": "Glutamatergic neuron",
"23": "Glutamatergic neuron",
"24": "Glutamatergic neuron",
"25": "Oligodendrocyte", 
"26": "Oligodendrocyte", 
"27": "Other neuron",
"28": "Glutamatergic neuron",
"29": "Astrocyte",
"30": "Cholinergic neuron",
"31": "Glutamatergic neuron",
"32": "Glutamatergic neuron",
"33": "GABAergic neuron",
"34": "Astrocyte",
"35": "Oligodendrocyte",
"36": "Glutamatergic neuron"
}


## UMAP with cell type labels
adata.obs['Cell type'] = adata.obs.leiden.map(Cell_type)
sc.pl.umap(adata, color = ['Cell type'])

## Set up dictionary for marker genes
marker_genes_dict = {
    'Glutamatergic neuron': ['Slc17a6', 'Slc17a7'],
    'GABAergic neuron': ['Gad1'],
    'Cholinergic neuron': ['Chat'],
    'Other neuron': ['Nos1', 'Cbln1', 'Cbln2', 'Bdnf'],
    'Microglia/Macrophage': ['Selplg', 'Slc15a3'],
    'Astrocyte': ['Aqp4', 'Aldh1l1', 'Cxcl14', 'Mlc1'],
    'Endothelial cell': ['Igf1r', 'Fn1', 'Slco1a4'],
    'Pericytes': ['Ace2', 'Adora2a', 'Rgs5'],
    'Oligodendrocytes': ['Gjc3', 'Opalin', 'Ermn', 'Sox8', 'Sgk1'],
}



##############################
# Build the spatial neighbors graph
sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)


# Compute centrality scores
sq.gr.centrality_scores(adata, cluster_key="Cell type")
sq.pl.centrality_scores(adata, cluster_key="Cell type", figsize=(16, 5))

# Co-occurrence probability
adata_subsample = sc.pp.subsample(adata, fraction=0.5, copy=True)

sq.gr.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
)
sq.pl.co_occurrence(
    adata_subsample,
    cluster_key="leiden",
    clusters="12",
    figsize=(10, 10),
)
sq.pl.spatial_scatter(
    adata_subsample,
    color="Cell type",
    shape=None,
    size=2,
)

# Neighborhood enrichment analysis
sq.gr.nhood_enrichment(adata, cluster_key="leiden")

ig, ax = plt.subplots(1, 2, figsize=(13, 7))
sq.pl.nhood_enrichment(
    adata,
    cluster_key="leiden",
    figsize=(8, 8),
    title="Neighborhood enrichment adata",
    ax=ax[0],
)
sq.pl.spatial_scatter(adata_subsample, color="leiden", shape=None, size=2, ax=ax[1])



####

# save figures as high-res
sc.set_figure_params(scanpy = True, dpi = 600, dpi_save = 600)
marker_genes_dict = {
    'Glutamatergic neuron': ['Slc17a6', 'Slc17a7'],
    'GABAergic neuron': ['Gad1'],
    'Cholinergic neuron': ['Chat'],
    'Other neuron': ['Nos1', 'Cbln1', 'Cbln2', 'Bdnf'],
    'Microglia/Macrophage': ['Selplg', 'Slc15a3'],
    'Astrocyte': ['Aqp4', 'Aldh1l1', 'Cxcl14', 'Mlc1'],
    'Endothelial cell': ['Igf1r', 'Fn1', 'Slco1a4'],
    'Pericytes': ['Ace2', 'Adora2a', 'Rgs5'],
    'Oligodendrocytes': ['Gjc3', 'Opalin', 'Ermn', 'Sox8', 'Sgk1'],
}
marker_genes = ['Slc17a6', 'Slc17a7', 'Gad1', 'Chat', 'Cbln1', 'Cbln2', 'Bdnf', 'Selplg', 'Slc15a3', 'Aldh1l1', 'Cxcl14', 'Mlc1', 'Igf1r', 'Fn1', 'Slco1a4', 'Adora2a', 'Rgs5', 'Gjc3', 'Opalin', 'Ermn', 'Sox8', 'Sgk1']


sc.pl.dotplot(adata, marker_genes, groupby='Cell type', save = 'cell_types.png')
sc.pl.umap(adata, color = ['leiden'], legend_loc = "on data", save = 'leiden_clusters.png')
sc.pl.umap(adata, color = ['Cell type'], save = 'cell_types.png')
sq.pl.spatial_scatter(adata, color="leiden", shape=None, axis_label = ['Spatial 1', 'Spatial 2'], legend_loc = None, dpi = 600, save = 'spatial_graph.png')
sq.pl.spatial_scatter(adata, color="leiden", legend_loc = None, axis_label = ['Spatial 1', 'Spatial 2'], save = 'spatial_graph.png')
sq.pl.nhood_enrichment(adata, cluster_key="leiden", mode = 'zscore', figsize=(8, 8), title="Neighborhood enrichment analysis",save = 'neighborhood_enrichment.png')

