#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 21:02:10 2024

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
import decoupler as dc

###############################################################################

# SETTINGS


## Matplotlib
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 500

###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/spatial'

adata = sc.read_h5ad(os.path.join(data_dir, 'adata-spatial.h5ad')) # immune cells
#adata = sc.read_h5ad(os.path.join(data_dir, 'bdata-spatial.h5ad')) # microglia

###############################################################################

# RUN DECOUPLER

progeny = dc.get_progeny(organism='mouse', top=500)
progeny


dc.run_mlm(
    mat=adata,
    net=progeny,
    source='source',
    target='target',
    weight='weight',
    verbose=True
)

adata.obsm['progeny_mlm_estimate'] = adata.obsm['mlm_estimate'].copy()
adata.obsm['progeny_mlm_pvals'] = adata.obsm['mlm_pvals'].copy()
adata

###############################################################################

# VISUALIZATION


acts = dc.get_acts(adata, obsm_key='mlm_estimate')
acts


sc.pl.umap(acts, color=['Trail', 'recluster'], cmap='RdBu_r', vcenter=0)
sc.pl.violin(acts, keys=['JAK-STAT'], groupby='recluster', rotation=90)


# PROGENY

sc.pl.matrixplot(acts, var_names=acts.var_names, groupby='recluster', dendrogram=False, standard_scale='var', colorbar_title='Z-scaled scores', cmap='viridis')

#sc.pl.matrixplot(acts, var_names=acts.var_names, groupby='cluster', dendrogram=False, standard_scale='var', colorbar_title='Z-scaled scores', cmap='viridis')


###############################################################################

# PATHWAY ENRICHMENT

msigdb = dc.get_resource('MSigDB')
msigdb
msigdb['collection'].unique()
print(msigdb['collection'].cat.categories)


# Filter by hallmark
msigdb = msigdb[msigdb['collection']=='hallmark']

# Remove duplicated entries
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
msigdb

dc.run_ora(
    mat=adata,
    net=msigdb,
    source='geneset',
    target='genesymbol',
    verbose=True,
    min_n = 2,
)

adata.obsm['ora_estimate']

###############################################################################

# VISUALIZATION

acts = dc.get_acts(adata, obsm_key='ora_estimate')

# We need to remove inf and set them to the maximum value observed
acts_v = acts.X.ravel()
max_e = np.nanmax(acts_v[np.isfinite(acts_v)])
acts.X[~np.isfinite(acts.X)] = max_e

acts

df = dc.rank_sources_groups(acts, groupby='recluster', reference='rest', method='wilcoxon')
df = pd.DataFrame(df)


n_markers = 3
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
source_markers

sc.pl.matrixplot(acts, source_markers, 'recluster', dendrogram=False, standard_scale='var')


###############################################################################


go_list = ['GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS',
           'GOBP_DEFENSE_RESPONSE',
           'GOBP_PROTEIN_MATURATION',
           'GOBP_PROTEOLYSIS']

sc.pl.matrixplot(acts, go_list, 'cluster', dendrogram=False, standard_scale='var', swap_axes = True, cmap = 'viridis')

sc.pl.matrixplot(acts, go_list, 'cluster', dendrogram=False, standard_scale='var', swap_axes = True, cmap = 'RdBu_r')


###############################################################################