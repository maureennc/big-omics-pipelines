#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 08:53:21 2024

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


###############################################################################

# SETTINGS


## Matplotlib
%matplotlib inline
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 8
plt.rcParams['figure.dpi'] = 500
plt.rcParams['figure.figsize'] = (3, 3)

###############################################################################

# IMPORT DATA

save_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))
                     
adata.obs.groupby('condition')['cell_type'].value_counts()

###############################################################################

# CELL COMPOSITION ANALYSIS - SAMPLE
df = pd.DataFrame(adata.obs)

df['Cluster'] = df['cluster'].astype('category')
df['Group'] = df['sample'].astype('category')

count_df = df.groupby(['Group', 'Cluster']).size().reset_index(name='Count')

plt.close()
g = sns.catplot(
    x='Cluster',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=2,
    palette='rocket',
    legend= True
)

g._legend.set_bbox_to_anchor((1.05, 1))  # (x, y) The position of the legend's bounding box
g._legend.set_title('sample') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()



df['Cell type'] = df['cell_type'].astype('category')
df['Group'] = df['sample'].astype('category')

count_df = df.groupby(['Group', 'Cell type']).size().reset_index(name='Count')

plt.close()
g = sns.catplot(
    x='Cell type',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=2,
    palette='rocket',
    legend= True
)

g._legend.set_bbox_to_anchor((1.05, 1))  # (x, y) The position of the legend's bounding box
g._legend.set_title('sample') 
plt.legend(title='Condition', loc='upper right', frameon=False)  # Move inside and remove border

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)

plt.show()

###############################################################################

# CELL TYPE COMPOSITION - GROUP

df = pd.DataFrame(adata.obs)

## Get value counts of cell types per condition using groupby and value_counts
cell_type_counts_by_condition = df.groupby('condition')['cell_type'].value_counts().reset_index(name='Count')

## Normalize counts by dividing by the number of samples per condition
sample_counts = {'Infected': 2, 'Naive': 2}

## Apply the normalization based on condition (divide by sample size for each condition)
cell_type_counts_by_condition['Normalized Count'] = cell_type_counts_by_condition.apply(
    lambda row: row['Count'] / sample_counts.get(row['condition'], 1), axis=1)

## Define the  order of cell types and conditions
cell_type_order = ['Excitatory neuron', 'Inhibitory neuron', 'Oligodendrocyte', 'Endothelial cell', 'Pericyte', 'Astrocyte', 'Microglia', 'Macrophage', 'CD4+ T cell', 'CD8+ T cell', 'Choroid plexus']

condition_order = ['Naive', 'Infected']  # Desired order for the condition hues

plt.figure(figsize=(3, 3))
sns.barplot(
    x='cell_type',
    y='Normalized Count',
    hue='condition',
    data=cell_type_counts_by_condition,
    palette='rocket',
    dodge=True,
    order=cell_type_order,  # Specify the order of the cell types on the x-axis
    hue_order=condition_order  # Specify the order of the conditions in the hue
)

plt.xticks(rotation=90)
plt.xlabel('')
plt.ylabel('Mean Cell Count')
plt.title('')
plt.legend(title='Condition', loc='upper right', frameon=False)
#plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.grid(False)
plt.show()

###############################################################################

# UMAP

plt.rcParams['font.size'] = 12

sc.pl.umap(adata, color = 'cluster', title = '')


## Control vs. infected umaps

condition_palette = ['#1B9E77', '#E7298A']
sc.pl.umap(adata, color = 'condition', title = '', palette = condition_palette)

infected_palette = ['#D73027', '#4575B4']  # red and blue
naive_palette = ['#66C2A5', '#FC8D62']  # green and orange

sc.pl.umap(adata, color='sample', groups=['E003-infected', 'E007-infected'], title='', palette=infected_palette)

sc.pl.umap(adata, color='sample', groups=['E008-naive', 'E009-naive'], title='', palette=naive_palette)

###############################################################################
