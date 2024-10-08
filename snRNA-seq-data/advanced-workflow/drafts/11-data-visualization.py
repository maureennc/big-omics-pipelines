#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 24 11:05:45 2024

@author: maureen
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import scanpy as sc


plt.rcParams['font.family'] = 'Arial'
sc.set_figure_params(scanpy = True, dpi = 300, dpi_save = 400, fontsize=15)

os.chdir('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/composition')

################################################################################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '3-tbi-annotated-hvg.h5ad'))
adata.X = adata.layers['log1p'].copy()

################################################################################################################################

# PREPARE DATA
save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/figures/summary/composition'

## Exclude artifacts
adata = adata[~adata.obs['cell_type'].isin(['Artifact', 'Unknown'])].copy()

## Update annotations
adata.obs['cell_type'] = adata.obs['cell_type'].replace('Unassigned', 'Unknown')
adata.obs['cell_class'] = adata.obs['cell_class'].replace('Unassigned', 'Unknown')
adata.obs['cluster'] = adata.obs['cluster'].replace('Unassigned', 'Unknown')

## Alphabetize annotations
adata.obs['cluster'] = pd.Categorical(adata.obs['cluster'], categories=sorted(adata.obs['cluster'].unique()), ordered=True) # commented 6/25/24
#adata.obs['cell_type'] = pd.Categorical(adata.obs['cell_type'], categories=sorted(adata.obs['cell_type'].unique()), ordered=True)


# Renumber leiden scores

value_counts = adata.obs['leiden_scVI'].value_counts()

cluster_mapping = {old_label: rank for rank, old_label in enumerate(value_counts.index)}
adata.obs['leiden'] = adata.obs['leiden_scVI'].map(cluster_mapping)
print(adata.obs['leiden'].value_counts())


################################################################################################################################

# CLUSTER DEFINING GENES - adata.obs.cluster

## Run wilcoxon rank-sum test
sc.tl.rank_genes_groups(adata, groupby='cluster', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 15, fontsize = 10, save = '_cluster-wilcoxon.pdf')

wilcoxon_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/cluster-gene-markers'

## Export gene rankings
markers_cluster = sc.get.rank_genes_groups_df(adata, None)
markers_cluster.to_csv(os.path.join(wilcoxon_dir, 'wilcoxon-cluster-markers-full.csv'), index = True) # export

## Extract top hits for each cluster
significant_markers = markers_cluster[(markers_cluster['pvals_adj'] < 0.05)].copy()
positive_markers = significant_markers[(significant_markers['logfoldchanges'] > 0)].copy()
top_genes_per_cluster = positive_markers.sort_values(by=['group', 'pvals_adj'], ascending=[True, True]).groupby('group').head(3)

## Plotting

sc.set_figure_params(scanpy = True, fontsize=16)

sc.tl.dendrogram(adata, groupby = 'cluster')
genes = top_genes_per_cluster['names'].tolist()

sc.pl.dotplot(adata, genes, groupby = 'cluster', dendrogram = False, standard_scale = 'var', save = 'cluster-marker-genes.pdf')
sc.pl.matrixplot(adata, genes, groupby = 'cluster', dendrogram = False, standard_scale = 'var', save = 'cluster-marker-genes.pdf')

################################################################################################################################

## CLUSTER DEFINING GENES - adata.obs.leiden


## Run wilcoxon rank-sum test
sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 15, fontsize = 10, save = '_leiden-cluster-wilcoxon.pdf')


## Export gene rankings
markers_leiden = sc.get.rank_genes_groups_df(adata, None)
markers_leiden.to_csv(os.path.join(wilcoxon_dir, 'wilcoxon-leiden-markers-full.csv'), index = True) # export

## Extract top hits for each leiden
significant_markers = markers_leiden[(markers_leiden['pvals_adj'] < 0.05)].copy()
positive_markers = significant_markers[(significant_markers['logfoldchanges'] > 0)].copy()
top_genes_per_leiden = positive_markers.sort_values(by=['group', 'pvals_adj'], ascending=[True, True]).groupby('group').head(3)

## Plotting

sc.set_figure_params(scanpy = True, fontsize=16)

sc.tl.dendrogram(adata, groupby = 'leiden')
genes = top_genes_per_leiden['names'].tolist()

sc.pl.dotplot(adata, genes, groupby = 'leiden', dendrogram = False, standard_scale = 'var', save = 'leiden-marker-genes.pdf')
sc.pl.matrixplot(adata, genes, groupby = 'leiden', dendrogram = False, standard_scale = 'var', save = 'leiden-marker-genes.pdf')

################################################################################################################################

# CREATE TABLE FOR LEIDEN LEGEND

csv_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/cluster-gene-markers'

leiden_legend = adata.obs[['leiden', 'cluster']].drop_duplicates()
leiden_legend = leiden_legend.sort_values(by = 'leiden')
leiden_legend.to_csv(os.path.join(csv_dir, 'leiden-legend.csv'), index = False)

################################################################################################################################

# CELL TYPE VALIDATION MARKERS

## Run wilcoxon rank-sum test
sc.tl.rank_genes_groups(adata, groupby='cell_type', method='wilcoxon', use_raw = True)
sc.pl.rank_genes_groups(adata, n_genes = 15, fontsize = 10, save = '_cell_type-wilcoxon.pdf')

## Extract top hits for each cell type
markers = sc.get.rank_genes_groups_df(adata, None)
positive_markers = markers[(markers['logfoldchanges'] > 0)]
top_genes_per_cell_type = positive_markers.sort_values(by=['group', 'pvals_adj'], ascending=[True, True]).groupby('group').head(3)

## Plotting
sc.tl.dendrogram(adata, groupby = 'cell_type')

### Wilcoxon result genes
genes = top_genes_per_cell_type['names'].unique().tolist()

sc.pl.dotplot(adata, genes, groupby = 'cell_type', dendrogram = False, save = 'cell_type-top3-genes.pdf')
sc.pl.dotplot(adata, genes, groupby = 'cell_type', dendrogram = True, save = 'cell_type-top3-genes-dendrogram.pdf')
sc.pl.matrixplot(adata, genes, groupby = 'cell_type', dendrogram = False, save = 'cell_type-top3-genes.pdf')

### Classical gene markers
classical_genes = ['Rbfox3', 'Nrgn', 'Slc17a7',
                   'Meg3', 'Gad1', 'Gad2', 
                   'Pdgfra', 'Cspg5', 'Sox10',
                   'Mbp', 'Mog', 'Plp1',
                   'Htr2c', 'Ttr', 'Gas6',
                   'Hexb', 'Cx3cr1', 'P2ry12',
                   'Slc1a2', 'Aqp4', 'Gfap',
                   'Vim', 'Lama1', 'Col1a1', 'Dcn']

sc.pl.dotplot(adata, classical_genes, groupby = 'cell_type', standard_scale = 'var', use_raw = True, dendrogram = True, save = '_cell-type-classical-markers-dendrogram.pdf')

################################################################################################################################

# COMPOSITIONAL ANALYSIS - CELL TYPE SUMMARY 

## Prepare data
bdata = adata[adata.obs['group'].isin(['A', 'B', 'C', 'E'])].copy()

## Focus on ipsilateral hemispheres
condition = {
    'A': "Sham + AAV-GFP",
    'B': "Sham + AAV-VEGFC",
    'C': "TBI + AAV-GFP",
    'E': "TBI + AAV-VEGFC"
}

bdata.obs['condition'] = bdata.obs.group.map(condition)

## Plotting
df = pd.DataFrame(bdata.obs)

df['Cell type'] = df['cell_type'].astype('category')
df['Group'] = df['condition'].astype('category')

count_df = df.groupby(['Group', 'Cell type']).size().reset_index(name='Count')

g = sns.catplot(
    x='Cell type',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=0.75,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((.97, 0.85))  # (x, y)
g._legend.set_title('Group') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.xlabel('')
plt.savefig(os.path.join(save_dir, "cell_type-counts-by-group.pdf"), format='pdf')
plt.show()


################################################################################################################################ 

# COMPOSITIONAL ANALYSIS - CLUSTER SUMMARY 

cluster_order = ['Astrocyte 1', 
         'Astrocyte 2', 
         'Astrocyte 3',
         'Excitatory neuron 1', 
         'Excitatory neuron 2', 
         'Excitatory neuron 3', 
         'Excitatory neuron 4', 
         'Excitatory neuron 5', 
         'Excitatory neuron 6', 
         'Excitatory neuron 7',
         'Excitatory neuron 8',
         'Excitatory neuron 9',
         'Excitatory neuron 10',
         'Excitatory neuron 11',
         'Excitatory neuron 12',
         'Excitatory neuron 13',
         'Excitatory neuron 14',
         'Excitatory neuron 15', 
         'Excitatory neuron 16',
         'Fibroblast', 
         'Inhibitory neuron 1',
         'Inhibitory neuron 2',
         'Inhibitory neuron 3', 
         'Inhibitory neuron 4',
         'Inhibitory neuron 5',
         'Inhibitory neuron 6',
         'Inhibitory neuron 7',
         'Inhibitory neuron 8',
         'Microglia',
         'OPC',
         'Oligodendrocyte',
         'Unknown']

df['Cluster'] = pd.Categorical(df['cluster'], categories=cluster_order, ordered=True)

count_df = df.groupby(['Group', 'Cluster']).size().reset_index(name='Count')

g = sns.catplot(
    x='Cluster',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=2,
    palette='viridis',
    legend= True
)

g._legend.set_bbox_to_anchor((0.95, 0.85))  # (x, y)
#g._legend.set_bbox_to_anchor((0.95, 0.85))  # (x, y)
g._legend.set_title('Group') 

plt.xticks(rotation=90)
plt.xlabel('')
plt.tight_layout()
plt.grid(False)
plt.savefig(os.path.join(save_dir, "cluster-counts-by-group.pdf"), format='pdf')
plt.show()

################################################################################################################################

# CHARACTERIZE UNKNOWN CLUSTER

unknown_cluster = markers_cluster[markers_cluster['group'] == 'Unknown'].copy()
significant_markers = unknown_cluster[unknown_cluster['pvals_adj'] < 0.05].copy()

positive_markers = significant_markers[(significant_markers['logfoldchanges'] > 0)].copy()
top_up_genes = positive_markers['names'].tolist()
top_up_genes = top_up_genes[:50].copy()

sc.pl.dotplot(adata, top_up_genes, groupby = 'cluster', use_raw = True, dendrogram = True, save = 'unknown-cluster-top-upreg.pdf')
sc.pl.matrixplot(adata, top_up_genes, groupby = 'cluster', use_raw = True, dendrogram = True, save = 'unknown-cluster-top-upreg.pdf')

################################################################################################################################

# UMAP

## Reset color mappings
print(bdata.uns.keys())

condition_colors = {
    "Sham + AAV-GFP": "#377eb8",
    "Sham + AAV-VEGFC": "#ff7f00",
    "TBI + AAV-GFP": "#4daf4a",
    "TBI + AAV-VEGFC": "#f781bf"
}

# Ensure that 'condition' is a categorical variable
bdata.obs['condition'] = bdata.obs['condition'].astype('category')

# Assign colors based on the categories in 'condition'
# The order of colors in the list should correspond to the order of categories
bdata.uns['condition_colors'] = [condition_colors[cat] for cat in bdata.obs['condition'].cat.categories]


sc.pl.umap(bdata, color = 'condition', groups = ['Sham + AAV-GFP', 'TBI + AAV-GFP'], save = "_treatment_condition_AC.pdf")
sc.pl.umap(bdata, color = 'condition', groups = ['Sham + AAV-GFP', 'Sham + AAV-VEGFC'], save = "_treatment_condition_AB.pdf")
sc.pl.umap(bdata, color = 'condition', groups = ['TBI + AAV-GFP', 'TBI + AAV-VEGFC'], save = "_treatment_condition_CE.pdf")

################################################################################################################################