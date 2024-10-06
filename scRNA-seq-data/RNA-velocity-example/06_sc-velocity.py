#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 10:36:07 2024

@author: maureen
"""

import os
import scanpy as sc
import anndata as ad
import scvelo as scv
from pathlib import Path
import pandas as pd
import numpy as np
import random
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns

###############################################################################

# SETTINGS

## Random seed
random.seed(0)
np.random.seed(0)

## Matplotlib
%matplotlib qt5
plt.rcParams['font.family'] = 'Arial'

## Scanpy
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400)
scv.settings.set_figure_params('scvelo', dpi = 300, dpi_save=400)


###############################################################################

# IMPORT DATA

## virtual environment = 'sc-seq'

data_dir = "/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/annotated-trained"
adata = sc.read_h5ad(os.path.join(data_dir, "pp-annotated-model_A.h5ad"))

## Restore raw counts to adata.X
if "log1p" in adata.layers:
    adata.X = adata.layers["log1p"].copy()
else:
    print("Log1p layer not found")

loom_dir = "/Users/maureen/Documents/experiments/cite-seq/data/loom"

naive_loom = scv.read(os.path.join(loom_dir, 'naive_cellranger_output_velocyto.loom'))
INF1_loom = scv.read(os.path.join(loom_dir, 'INF1_cellranger_output_velocyto.loom'))
INF2_loom = scv.read(os.path.join(loom_dir, 'INF2_cellranger_output_velocyto.loom'))

naive_loom.var_names_make_unique()
INF1_loom.var_names_make_unique()
INF2_loom.var_names_make_unique()

###############################################################################

# EDIT BARCODES FOR MAPPING

## Naive
original_index = naive_loom.obs.index
processed_index = [barcode.replace('naive_cellranger_output_velocyto:', '').replace('x', '') for barcode in original_index]
naive_loom.obs.index = processed_index

## INF1
original_index = INF1_loom.obs.index
processed_index = [barcode.replace('INF1_cellranger_output_velocyto:', '').replace('x', '') for barcode in original_index]
INF1_loom.obs.index = processed_index

## INF2
original_index = INF2_loom.obs.index
processed_index = [barcode.replace('INF2_cellranger_output_velocyto:', '').replace('x', '') for barcode in original_index]
INF2_loom.obs.index = processed_index

###############################################################################

# CONCATENATE

ldata = naive_loom.concatenate(INF1_loom, INF2_loom, batch_key='sample', batch_categories=['naive', 'INF1', 'INF2'] )
ldata.obs

###############################################################################

# SUBSET ON INTERSECTION OF ADATA AND LOOM

## Merge    
vdata = scv.utils.merge(adata, ldata).copy()

## Subset out microglia
vdata = vdata[vdata.obs['cell_class'] == 'Microglia'].copy()

sc.pp.neighbors(vdata, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(vdata, min_dist = 0.3)
sc.tl.leiden(vdata, key_added='leiden_scVI', resolution=0.3)
sc.pl.umap(vdata, color = 'leiden_scVI')


cluster= { 
"0": "Mg-0",
"1": "Mg-3",
"2": "Mg-2",
"3": "Mg-1"
}

vdata.obs['cluster'] = vdata.obs.leiden_scVI.map(cluster)

order = ['Mg-0', 'Mg-1', 'Mg-2', 'Mg-3']
vdata.obs['cluster'] = pd.Categorical(vdata.obs['cluster'], categories=order, ordered=True)

sc.pl.umap(vdata, color = ['cluster'], legend_loc = 'on data', legend_fontsize=12, legend_fontoutline=2,)


###############################################################################

# REGRESS OUT EFFECT OF TOTAL COUNTS
sc.pp.regress_out(vdata, ['total_counts'])

###############################################################################

# TRAJECTORY INFERENCE ANALYSIS

scv.pl.proportions(vdata, groupby = 'cluster')

scv.pp.filter_and_normalize(vdata, log = False)

scv.pp.moments(vdata, n_pcs=20, n_neighbors=20)
scv.tl.recover_dynamics(vdata)

scv.tl.velocity(vdata, mode='dynamical')
scv.tl.velocity_graph(vdata)
scv.pl.velocity_embedding_stream(vdata, basis='umap', color = 'cluster')

## Latent time incorporates more information
scv.tl.latent_time(vdata)
scv.pl.scatter(vdata, color='latent_time', color_map='gnuplot', size=50,frameon = True, dpi = 300, rescale_color = (0, 1))
scv.pl.scatter(vdata, color='latent_time', color_map='gnuplot', size=20,frameon = True, dpi = 300, linewidth = 0.5, colorbar = False)

#scv.tl.velocity_pseudotime(vdata)
#scv.pl.scatter(vdata, color='velocity_pseudotime', cmap='gnuplot')


#scv.tl.terminal_states(vdata)
#scv.pl.scatter(vdata, color=['root_cells', 'end_points'])

scv.pl.velocity_embedding(vdata, basis='umap', arrow_length=3, arrow_size=2, dpi=120, color='latent_time')

## PAGA Embedding
vdata.uns['neighbors']['distances'] = vdata.obsp['distances']
vdata.uns['neighbors']['connectivities'] = vdata.obsp['connectivities']

scv.tl.paga(vdata, groups='cluster')
df = scv.get_df(vdata, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(vdata, basis='umap', size=50, alpha=.1, min_edge_width=2, node_size_scale=1.5)


scv.pl.velocity(vdata, ['Il1a', 'Nlrp1b'])

###############################################################################

# DATA VIZ

scv.pp.filter_and_normalize(vdata, log=False)
scv.pp.moments(vdata, n_pcs=20, n_neighbors=20)
scv.tl.recover_dynamics(vdata)
scv.tl.velocity(vdata, mode='dynamical')
scv.tl.velocity_graph(vdata)
scv.tl.latent_time(vdata)

# Visualize UMAP embedding with velocity streams and color by latent time
scv.pl.velocity_embedding_stream(
    vdata,
    basis='umap',  # Ensure you've computed UMAP previously
    color='latent_time',
    palette='gnuplot')


###############################################################################