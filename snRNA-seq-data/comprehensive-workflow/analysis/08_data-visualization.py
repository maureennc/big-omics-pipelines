#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  7 17:05:23 2024

@author: maureen

Figure 4 & Supplemental Figure 4

This script visualizes:
    - UMAPs
    - Differential expression summary heatmaps
    - Cell type composition summary
    - Heatmaps for genes associated with neuronal health and pathology
    - Violin plots for microglial and astrocytic genes of interest

## Additional visualizations in `collaborations` repo: `11-data-visualization.py`, under project folder
## Just including what we are submitting for publication here

"""

import os
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from functools import reduce
from functools import reduce

sc.set_figure_params(scanpy = True)


plt.rcParams['figure.dpi'] = 600
plt.rcParams['font.family'] = 'Arial'


###############################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '6-tbi-annotated-full.h5ad'))

bdata = sc.read_h5ad(os.path.join('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad/5-neuron-recluster-full.h5ad'))

###############################################################################

# FIGURE 4

###############################################################################

# UMAP

## Change directory for scanpy save figures to work
os.chdir('/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/submission_01')

## With legend
sc.pl.umap(adata, color = ['cell_type'], legend_fontsize = 14, title = 'Cell types', save = '_cell_type.pdf')
sc.pl.umap(bdata, color = ['neuron_cluster'], legend_fontsize = 10, title = 'Neuronal clusters', save = '_neuron_cluster.pdf')

## No legend
sc.pl.umap(adata, color=['cell_type'], legend_loc=None, title='Cell types', save='_cell_type_no_legend.pdf')

sc.pl.umap(bdata, color=['neuron_cluster'], legend_loc=None, title='Neuronal clusters', save='_neuron_cluster_no_legend.pdf')

###############################################################################

# CELL TYPE DIFFERENTIAL EXPRESSION SUMMARY HEATMAP

## Import DE results

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/3/results/de/cell_type"

final_results =  pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type.csv'))
filtered_results =  pd.read_csv(os.path.join(save_dir, 'wald-test-cell_type-filtered.csv'))


## Run differential expression

comparisons = ['AC', 'CE', 'AB']

order = ['Astrocyte', 'Excitatory neuron', 'Inhibitory neuron',  'Microglia', 'OPC', 'Oligodendrocyte', 'Fibroblast', 'Choroid-plexus epithelial' ]

## Create a reference DataFrame for cell types
cell_type_reference = pd.DataFrame(order, columns=['cell_type'])

data_list_up = []
data_list_down = []

mean_expression_threshold = 0.1
## Loop over each comparison
for comp in comparisons:
    df_comp = filtered_results[(filtered_results['comparison'] == comp)]
    df_comp = df_comp[df_comp['mean'] > mean_expression_threshold]

    # Upregulated genes
    df_up = df_comp[(df_comp['qval'] < 0.05) & (df_comp['log2fc'] > 0.5)]
    deg_counts_up = df_up.groupby('cell_type').size().reset_index(name=comp)  # Retain original comparison name
    deg_counts_up = pd.merge(cell_type_reference, deg_counts_up, on='cell_type', how='left').fillna(0)
    data_list_up.append(deg_counts_up)

    # Downregulated genes
    df_down = df_comp[(df_comp['qval'] < 0.05) & (df_comp['log2fc'] < -0.5)]
    deg_counts_down = df_down.groupby('cell_type').size().reset_index(name=comp)  # Retain original comparison name
    deg_counts_down = pd.merge(cell_type_reference, deg_counts_down, on='cell_type', how='left').fillna(0)
    data_list_down.append(deg_counts_down)

## Combine up and down dfs
combined_df_up = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list_up)
combined_df_down = reduce(lambda left, right: pd.merge(left, right, on='cell_type', how='outer'), data_list_down)
combined_df_up.fillna(0, inplace=True)
combined_df_down.fillna(0, inplace=True)

## Calculate vmin and vmax for 
vmin = min(combined_df_up.drop(columns='cell_type').min().min(), combined_df_down.drop(columns='cell_type').min().min())
vmax = max(combined_df_up.drop(columns='cell_type').max().max(), combined_df_down.drop(columns='cell_type').max().max())

# Reorder by cell type
combined_df_up['cell_type'] = pd.Categorical(combined_df_up['cell_type'], categories=order, ordered=True)
combined_df_up = combined_df_up.sort_values('cell_type')

combined_df_down['cell_type'] = pd.Categorical(combined_df_down['cell_type'], categories=order, ordered=True)
combined_df_down = combined_df_down.sort_values('cell_type')


save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/submission_01/figures'

## Upregulated
plt.figure(figsize=(6, 3))
sns.heatmap(combined_df_up.set_index('cell_type'), annot=True, cmap='viridis', fmt="g", vmin=vmin, vmax=vmax)
plt.title('Upregulated Genes')
plt.ylabel(' ')
plt.xlabel('')  # No label for the x-axis
plt.xticks(ticks=range(len(comparisons)), labels=comparisons)  # Set custom x-axis labels
plt.grid(False)
plt.tight_layout()
plt.xticks(ticks=np.arange(len(comparisons)) + 0.5, labels=comparisons, rotation=0)
plt.savefig(os.path.join(save_dir, 'upregulated-genes-summary.pdf'))
plt.show()

## Downregulated
plt.figure(figsize=(6, 3))
sns.heatmap(combined_df_down.set_index('cell_type'), annot=True, cmap='viridis', fmt="g", vmin=vmin, vmax=vmax)
plt.title('Downregulated Genes')
plt.ylabel(' ')
plt.xlabel('')  # No label for the x-axis
plt.xticks(ticks=range(len(comparisons)), labels=comparisons)  # Set custom x-axis labels
plt.grid(False)
plt.tight_layout()
plt.xticks(ticks=np.arange(len(comparisons)) + 0.5, labels=comparisons, rotation=0)
plt.savefig(os.path.join(save_dir, 'downregulated-genes-summary.pdf'))
plt.show()

###############################################################################

# COMPOSITIONAL ANALYSIS - CELL TYPE SUMMARY 

## Prepare data
adata = adata[adata.obs['group'].isin(['A', 'B', 'C', 'E'])].copy()

## Focus on ipsilateral hemispheres
condition = {
    'A': r'Sham + AAV$^{\mathrm{GFP}}$',
    'B': r'Sham + AAV$^{\mathrm{VEGFC}}$',
    'C': r'TBI + AAV$^{\mathrm{GFP}}$',
    'E': r'TBI + AAV$^{\mathrm{VEGFC}}$'
}

adata.obs['condition'] = adata.obs.group.map(condition)


### Add superscript for GFP/VEGFC annotation; makes it italicized
#condition = {
#    'A': r'Sham + AAV$^{GFP}$',
#    'B': r'Sham + AAV$^{VEGFC}$',
#    'C': r'TBI + AAV$^{GFP}$',
#    'E': r'TBI + AAV$^{VEGFC}$'
#}


## Plotting
save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/submission_01/figures'

df = pd.DataFrame(adata.obs)

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


###############################################################################

# NEURON SUBSET - CLUSTER COMPOSITION

groups = ['A', 'B', 'C', 'E']
cdata = bdata[bdata.obs['group'].isin(groups)].copy()

## Rename group annotations for this subset of bdata (neurons only)

### Not italicized
condition = {
    'A': r'Sham + AAV$^{\mathrm{GFP}}$',
    'B': r'Sham + AAV$^{\mathrm{VEGFC}}$',
    'C': r'TBI + AAV$^{\mathrm{GFP}}$',
    'E': r'TBI + AAV$^{\mathrm{VEGFC}}$'
}

cdata.obs['condition'] = cdata.obs['group'].map(condition)

   
## Bar chart  

plt.rcParams['font.size'] = 12
         
cluster_counts = cdata.obs.groupby(['condition', 'neuron_cluster']).size().reset_index(name='Count')
print(cluster_counts)

## Calculate the total number of neurons in each group
total_neurons_per_group = cdata.obs.groupby('condition').size().reset_index(name='Total_Neurons_group')
print(total_neurons_per_group)

## Merge the total neuron counts into the cluster counts DataFrame
cluster_counts = pd.merge(cluster_counts, total_neurons_per_group, on='condition')
cluster_counts

## Calculate the frequency of each neuron cluster within each group
cluster_counts['Frequency'] = cluster_counts['Count'] / cluster_counts['Total_Neurons_group']
cluster_counts
## Plot the cluster frequencies using seaborn
plt.figure(figsize=(3, 5))
g = sns.catplot(
    x='neuron_cluster', 
    y='Frequency',  # Use the calculated frequency
    hue='condition', 
    data=cluster_counts, 
    kind='bar', 
    height=3, 
    aspect=2, 
    palette='viridis'
)

# Customize the plot
g.set_axis_labels(" ", "% Neuronal Population\nper Group")  # Custom x and y labels
g._legend.set_title('Group')
g._legend.set_bbox_to_anchor((0.95, .9))  # Adjust the legend position
plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.savefig("./figures/cluster_frequencies.pdf", format='pdf', bbox_inches='tight') 
plt.show()


###############################################################################

## PATHOLOGY-ASSOCIATED GENES HEATMAP (FROM HOLTZMANN PAPER)

plt.rcParams['font.size'] = 16

## With Zbtb20
genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Brinp3', 'Mef2c', 'Zbtb20']

### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups-Zbtb20.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups-Zbtb20.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups-Zbtb20.pdf')

sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups-Zbtb20.pdf')


## Without Zbtb20
genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Brinp3', 'Mef2c']

### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-six-groups.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups.pdf')

sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'pathology-genes-four-groups.pdf')


###############################################################################

# NEURONAL HEALTH-ASSOCIATED GENES HEATMAP (HOLTZMANN PAPER)

genes = ['Prox1', 'Synpr', 'C1ql2', 'C1ql3', 'Camk2a', 'Camk2b', 'Tmem108', 'Ppfia2', 'Rfx3', 'Lrrtm4', 'Btbd9', 'Cntnap5a', 'Erc2']

plt.rcParams['font.size'] = 16


### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-six-groups.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-six-groups.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-four-groups.pdf')

sc.pl.dotplot(cdata, var_names = genes, groupby = 'neuron_cluster', standard_scale = 'var', dendrogram = True, swap_axes = True, save = 'neuronal-health-genes-four-groups.pdf',  cmap = 'viridis')


###############################################################################

# PATHOLOGY-ASSOCIATED GENES BY TREATMENT GROUP

genes = ['Arpp21', 'R3hdm1', 'Rorb', 'Cux1', 'Cux2', 'Brinp3', 'Mef2c', 'Zbtb20']

plt.rcParams['font.size'] = 16


### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var', dendrogram = True,  save = 'pathology-genes-six-groups-condition.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var',  save = 'pathology-genes-six-groups-condition.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'condition', standard_scale = 'var',  save = 'pathology-genes-four-groups-condition.pdf')

sc.pl.dotplot(cdata, var_names = genes, groupby = 'condition', standard_scale = 'var', save = 'pathology-genes-four-groups-condition.pdf',  cmap = 'viridis')


###############################################################################

# HEALTH-ASSOCIATED GENES BY TREATMENT GROUP

genes = ['Prox1', 'Synpr', 'C1ql2', 'C1ql3', 'Camk2a', 'Camk2b', 'Tmem108', 'Ppfia2', 'Rfx3', 'Lrrtm4', 'Btbd9', 'Cntnap5a', 'Erc2']

plt.rcParams['font.size'] = 16


### Six groups
sc.pl.matrixplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var', dendrogram = True,  save = 'neuronal-health-genes-six-groups-condition.pdf')

sc.pl.dotplot(bdata, var_names = genes, groupby = 'condition', standard_scale = 'var',  save = 'neuronal-health-genes-six-groups-condition.pdf', cmap = 'viridis')

### Four groups
sc.pl.matrixplot(cdata, var_names = genes, groupby = 'condition', standard_scale = 'var',  save = 'neuronal-health-genes-four-groups-condition.pdf')

sc.pl.dotplot(cdata, var_names = genes, groupby = 'condition', standard_scale = 'var', save = 'neuronal-health-genes-four-groups-condition.pdf',  cmap = 'viridis')

###############################################################################
###############################################################################

# SUPPLEMENTAL FIGURE 4

###############################################################################

# UMAPs

sc.pl.umap(adata, color = 'cluster', title = 'Cluster')

###############################################################################

# DOTPLOT ALL CLUSTERS

sc.tl.rank_genes_groups(adata, groupby = 'cluster', method = 'wilcoxon')
cluster_markers = sc.get.rank_genes_groups_df(adata, group = None)
cluster_markers

## Gene list - for plotting top 3 genes for each cluster
gene_list = cluster_markers.groupby('group').head(3)
genes = gene_list['names'].tolist()

sc.pl.dotplot(adata, var_names = genes, groupby = 'cluster')

###############################################################################

# VIOLIN PLOTS

plt.rcdefaults()
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 14
plt.rcParams['figure.dpi'] = 600

## Astrocytes

# astrocyte_genes = ['Gfap', 'Apoe', 'S100b']

astrocyte = adata[adata.obs['cell_type'] == 'Astrocyte'].copy()
groups = ['A', 'B', 'C', 'E']

astrocyte = astrocyte[astrocyte.obs['group'].isin(groups)]

condition = {
    'A': r'Sham + AAV$^{\mathrm{GFP}}$',
    'B': r'Sham + AAV$^{\mathrm{VEGFC}}$',
    'C': r'TBI + AAV$^{\mathrm{GFP}}$',
    'E': r'TBI + AAV$^{\mathrm{VEGFC}}$'
}

astrocyte.obs['condition'] = astrocyte.obs['group'].map(condition)

## Group
plt.rcParams['figure.figsize'] = [4, 3] 

sc.pl.violin(astrocyte, keys = ['Gfap'], groupby = 'group', rotation = 45, use_raw = True, save = '_astrocyte-group-Gfap.pdf')
sc.pl.violin(astrocyte, keys = ['Apoe'], groupby = 'group', rotation = 45, use_raw = True, save = '_astrocyte-group-Apoe.pdf')
sc.pl.violin(astrocyte, keys = ['S100b'], groupby = 'group', rotation = 45, use_raw = True, save = '_astrocyte-group-S100b.pdf')

## Condition
plt.rcParams['figure.figsize'] = [4, 2] 

sc.pl.violin(astrocyte, keys = ['Gfap'], groupby = 'condition', rotation = 45, use_raw = True, save = '_astrocyte-condition-Gfap.pdf')
sc.pl.violin(astrocyte, keys = ['Apoe'], groupby = 'condition', rotation = 45, use_raw = True, save = '_astrocyte-condition-Apoe.pdf')
sc.pl.violin(astrocyte, keys = ['S100b'], groupby = 'condition', rotation = 45, use_raw = True, save = '_astrocyte-condition-S100b.pdf')



## Microglia
# microglia_genes = ['Cd68', 'H2-Ab1', 'P2ry12']

microglia = adata[adata.obs['cell_type'] == 'Microglia'].copy()
groups = ['A', 'B', 'C', 'E']

microglia = microglia[microglia.obs['group'].isin(groups)]

condition = {
    'A': r'Sham + AAV$^{\mathrm{GFP}}$',
    'B': r'Sham + AAV$^{\mathrm{VEGFC}}$',
    'C': r'TBI + AAV$^{\mathrm{GFP}}$',
    'E': r'TBI + AAV$^{\mathrm{VEGFC}}$'
}

microglia.obs['condition'] = microglia.obs['group'].map(condition)

## Group
plt.rcParams['figure.figsize'] = [4, 3] 

sc.pl.violin(microglia, keys = ['Cd68'], groupby = 'group', rotation = 45, use_raw = True, save = '_microglia-group-Cd68.pdf')
sc.pl.violin(microglia, keys = ['H2-Ab1'], groupby = 'group', rotation = 45, use_raw = True, save = '_microglia-group-H2-Ab1.pdf')
sc.pl.violin(microglia, keys = ['P2ry12'], groupby = 'group', rotation = 45, use_raw = True, save = '_microglia-group-P2ry12.pdf')

## Condition
plt.rcParams['figure.figsize'] = [4, 2] 

sc.pl.violin(microglia, keys = ['Cd68'], groupby = 'condition', rotation = 45, use_raw = True, save = '_microglia-condition-Cd68.pdf')
sc.pl.violin(microglia, keys = ['H2-Ab1'], groupby = 'condition', rotation = 45, use_raw = True, save = '_microglia-condition-H2-Ab1.pdf')
sc.pl.violin(microglia, keys = ['P2ry12'], groupby = 'condition', rotation = 45, use_raw = True, save = '_microglia-condition-P2ry12.pdf')

###############################################################################
