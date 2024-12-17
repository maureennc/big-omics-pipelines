#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 16 13:53:54 2024

@author: maureen
"""


import os
import scanpy as sc
import decoupler as dc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 16

plt.rcParams['figure.dpi'] = 300


###############################################################################

# SETTINGS

#sc.set_figure_params(scanpy = True, dpi = 100, dpi_save = 200, fontsize = 14, figsize = None)

###############################################################################

# IMPORT DATA

## adata

data_dir = '/Users/maureen/Documents/experiments/cite-seq/analysis/2/python/h5ad/all-cells/pre-processed'

adata = sc.read_h5ad(os.path.join(data_dir, "pp-concat-full-metadata.h5ad"))

###############################################################################

# Prepare data

## Subset out microglia
adata = adata[adata.obs['cell_class'] == 'Microglia'].copy()
print(adata.X[:5])  # Print the first five entries

adata.X = adata.layers['counts'].copy()
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)

sc.pp.regress_out(adata, keys='total_counts')
sc.pp.scale(adata)

sc.pp.neighbors(adata, use_rep = 'X_scVI', random_state = 0) # Use latent representation to build neighbors graph
sc.tl.umap(adata, min_dist = 0.3)
sc.tl.leiden(adata, key_added='leiden_scVI', resolution=0.3)
sc.pl.umap(adata, color = 'leiden_scVI')


cluster= { 
"0": "Mg-0",
"1": "Mg-3",
"2": "Mg-2",
"3": "Mg-1"
}

adata.obs['cluster'] = adata.obs.leiden_scVI.map(cluster)

order = ['Mg-0', 'Mg-1', 'Mg-2', 'Mg-3']

adata.obs['cluster'] = pd.Categorical(adata.obs['cluster'], categories=order, ordered=True)
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


sc.pl.umap(acts, color=['Trail', 'cluster'], cmap='RdBu_r', vcenter=0)
sc.pl.violin(acts, keys=['MAPK'], groupby='cluster', rotation=90)


# PROGENY

sc.pl.matrixplot(acts, var_names=acts.var_names, groupby='cluster', dendrogram=False, standard_scale='var', colorbar_title='Z-scaled scores', cmap='viridis')

#sc.pl.matrixplot(acts, var_names=acts.var_names, groupby='cluster', dendrogram=False, standard_scale='var', colorbar_title='Z-scaled scores', cmap='viridis')



genes_list = [
    'Apoe', 'Itgax', 'Axl', 'Ccl2', 'Tlr2','Cybb', 'Csf1', 'Lag3',
   'Gpx3', 'Ccrl2', 'Cxcl10', 'Cxcl16', 'Cxcr4',
    'Gpnmb', 'Lgals3', 
    'Spp1', 'Msr1', 'Arg1', 'Cfp', 'Alcam', 'Gas7', 'Siglec1' #'Ifi202b'
]


sc.pl.matrixplot(adata, var_names=genes_list, groupby='cluster', dendrogram=False, standard_scale='var', colorbar_title='Z-scaled scores' )

sc.pl.heatmap(adata, var_names=genes_list, groupby='group', dendrogram=False)

parasite_genes = [
    'Gbp2b', 'Gbp2', 'Gbp3', 'Gbp4', 'Gbp6', 'Gbp8', 'Gbp9', 'Gbp10', 
    'Irgm2', 'Iigp1', 'Il12b', 'Myd88', 'Cd40', 'Slc11a1', 'Batf2', 
    'Irf8', 'Tlr12', 'Ier3', 'Il4ra', 'Pf4', 'Arg1', 'Enpp1', 'Il6', 'Irf4', 
]

sc.pl.matrixplot(adata, var_names=parasite_genes, groupby='cluster', dendrogram=False, standard_scale='var', colorbar_title='Z-scaled scores', cmap='viridis')



ginhoux_genes = ['Dkk2', 'Fabp5', 'Gpnmb', 'Igf1', 'Itgax', 'Mamdc2', 'Spp1', 'Gm1673']
sc.pl.matrixplot(adata, var_names=ginhoux_genes, groupby='cluster', dendrogram=False, standard_scale='var')
###############################################################################

# PATHWAY ENRICHMENT

msigdb = dc.get_resource('MSigDB')
msigdb
msigdb['collection'].unique()
print(msigdb['collection'].cat.categories)


# Filter by hallmark
msigdb = msigdb[msigdb['collection']=='go_biological_process']

# Remove duplicated entries
msigdb = msigdb[~msigdb.duplicated(['geneset', 'genesymbol'])]
msigdb

dc.run_ora(
    mat=adata,
    net=msigdb,
    source='geneset',
    target='genesymbol',
    verbose=True
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

df = dc.rank_sources_groups(acts, groupby='cluster', reference='rest', method='wilcoxon')
df = pd.DataFrame(df)


n_markers = 5
source_markers = df.groupby('group').head(n_markers).groupby('group')['names'].apply(lambda x: list(x)).to_dict()
source_markers

sc.pl.matrixplot(acts, source_markers, 'cluster', dendrogram=True, standard_scale='var')


###############################################################################


go_list = ['GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS',
           'GOBP_DEFENSE_RESPONSE',
           'GOBP_PROTEIN_MATURATION',
           'GOBP_PROTEOLYSIS']

sc.pl.matrixplot(acts, go_list, 'cluster', dendrogram=False, standard_scale='var', swap_axes = True, cmap = 'viridis')

sc.pl.matrixplot(acts, go_list, 'cluster', dendrogram=False, standard_scale='var', swap_axes = True, cmap = 'RdBu_r')


###############################################################################

# Differential expression

sc.tl.rank_genes_groups(adata, groupby = 'cluster')
sc.pl.rank_genes_groups(adata, groupby = 'cluster')

de = sc.get.rank_genes_groups_df(adata, group = None)

save_dir = '/Users/maureen/Desktop/allen-institute/tables'
de.to_csv(os.path.join(save_dir, 'sc-seq-microglia-clusters-markers.csv'), index=False)

# Ginhoux genes
#genes = ['Dkk2', 'Fabp5', 'Gpnmb', 'Igf1', 'Itgax', 'Mamdc2', 'Spp1', 'Gm1673']
#sc.pl.matrixplot(adata, var_names = genes, groupby = 'cluster', standard_scale = 'var')

###############################################################################

# PIE CHARTS

plt.rcParams['font.size'] = 16

# Define a more appealing color mapping for your clusters
cluster_colors = {
    'Mg-0': '#1f77b4',  # Color for Mg-0 (blue)
    'Mg-1': '#ff7f0e',  # Color for Mg-1 (orange)
    'Mg-2': '#2ca02c',  # Color for Mg-2 (green)
    'Mg-3': '#d62728'   # Color for Mg-3 (red)
}

# Get counts for each cluster in naive and infected groups
naive_cluster_counts = naive_cells.obs['cluster'].value_counts()
infected_cluster_counts = infected_cells.obs['cluster'].value_counts()

# Ensure the colors are in the same order as the counts
naive_colors = [cluster_colors[cluster] for cluster in naive_cluster_counts.index]
infected_colors = [cluster_colors[cluster] for cluster in infected_cluster_counts.index]

# Function to format the label based on value
def func(pct, allvalues, threshold=3):
    return f'{pct:.1f}%' if pct >= threshold else ''  # Show percentage only if it is greater than or equal to the threshold

# Create pie chart for naive group
plt.figure(figsize=(10, 5))

plt.subplot(1, 2, 1)  # 1 row, 2 columns, 1st subplot
wedges, texts, autotexts = plt.pie(
    naive_cluster_counts, 
    autopct=lambda pct: func(pct, naive_cluster_counts, threshold=5),  # Use the custom function with threshold
    startangle=90,
    colors=naive_colors,  # Use the specific colors for naive
    textprops=dict(color="black"),  # Set font properties for wedge labels
    explode=(0.1,) * len(naive_cluster_counts),  # Slightly explode slices for visibility
    pctdistance=0.6  # Move the percentage labels back toward the center
)

plt.title(' ')
plt.axis('equal')  # Equal aspect ratio ensures that pie chart is circular.

# Create pie chart for infected group
plt.subplot(1, 2, 2)  # 1 row, 2 columns, 2nd subplot
wedges, texts, autotexts = plt.pie(
    infected_cluster_counts, 
    autopct=lambda pct: func(pct, infected_cluster_counts, threshold=5),  # Use the custom function with threshold
    startangle=90,
    colors=infected_colors,  # Use the specific colors for infected
    textprops=dict(color="black"),  # Set font properties for wedge labels
    explode=(0.1,) * len(infected_cluster_counts),  # Slightly explode slices for visibility
    pctdistance=0.6  # Move the percentage labels back toward the center
)

plt.title(' ')
plt.axis('equal')  # Equal aspect ratio ensures that pie chart is circular.

# Create a legend for both pie charts with increased font size
handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10) for color in naive_colors]
plt.figlegend(handles, naive_cluster_counts.index, title="Clusters", loc="center right", bbox_to_anchor=(1.1, 0.5), fontsize='large')  # Increase font size

# Show the plots
plt.tight_layout()
plt.show()


###############################################################################

# visualizations

sc.pl.violin(adata, keys = ['total_counts', 'Trem2'], groupby = 'cluster')

sc.pl.violin(adata, keys = ['Apoe', 'Trem2', 'Cd33'], groupby = 'cluster', palette = 'tab10')

#sc.pl.violin(adata, keys = ['Apoe', 'Trem2', 'Tyrobp', 'Cd33'], groupby = 'cluster', palette = 'tab10')

sc.pl.violin(adata, keys = ['total_counts', 'n_genes_by_counts', 'pct_counts_ribosomal'], groupby = 'cluster', palette = 'tab10')

sc.pl.violin(adata, keys = ['Gapdh', 'Jun', 'Arc'], groupby = 'cluster')
