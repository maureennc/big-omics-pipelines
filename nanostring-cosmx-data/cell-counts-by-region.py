#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 12:06:05 2024

@author: maureen
"""

import os
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scanpy as sc
import random
import torch

###############################################################################

# SETTINGS

## Random seed
random.seed(0)
torch.manual_seed(0)
np.random.seed(0)

###############################################################################

# IMPORT DATA 
data_dir = "/Users/maureen/Documents/projects/lukens-lab/nick/syk-aging-cosmx/analysis/3/h5ad/"

adata = sc.read_h5ad(os.path.join(data_dir, '2-syk-integrated-annotated.h5ad'))

###############################################################################

# CELL COMPOSITION ANALYSIS

## Cell type
df = pd.DataFrame(adata.obs)

df['Cell type'] = df['cell_type'].astype('category')
df['Group'] = df['region'].astype('category')

count_df = df.groupby(['Group', 'Cell type']).size().reset_index(name='Count')

plt.close()
g = sns.catplot(
    x='Cell type',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=.75,
    palette='rocket',
    legend= True
)

g._legend.set_bbox_to_anchor((1.05, 1))
g._legend.set_title('Region') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()


## Cluster
df = pd.DataFrame(adata.obs)

df['Cluster'] = df['cluster'].astype('category')
df['Group'] = df['region'].astype('category')

count_df = df.groupby(['Group', 'Cluster']).size().reset_index(name='Count')

plt.close()
g = sns.catplot(
    x='Cluster',
    y='Count',
    hue='Group',  
    data=count_df,
    kind='bar',
    height=5,
    aspect=.75,
    palette='rocket',
    legend= True
)

g._legend.set_bbox_to_anchor((1.05, 1))
g._legend.set_title('Dataset') 

plt.xticks(rotation=90)
plt.tight_layout()
plt.grid(False)
plt.show()

###############################################################################

# CELL COUNTS BY REGION
groups = ['Oligodendrocyte', 'Microglia']

hippocampus_data = adata[adata.obs['region'] == 'Hippocampus'].copy()
hippocampus_data = pd.DataFrame(hippocampus_data.obs)
hippocampus_data = hippocampus_data[hippocampus_data['cell_type'].isin(groups)]
hippocampus_data['cell_type'] = hippocampus_data['cell_type'].cat.remove_unused_categories()

callosum_data = adata[adata.obs['region'] == 'Corpus callosum'].copy()
callosum_data = pd.DataFrame(callosum_data.obs)
callosum_data = callosum_data[callosum_data['cell_type'].isin(groups)]
callosum_data['cell_type'] = callosum_data['cell_type'].cat.remove_unused_categories()

## Hippocampus

hippo_counts = hippocampus_data.groupby(['cell_type', 'condition'], observed=True).size().reset_index(name='counts').copy()

plt.figure(figsize=(7, 5), dpi = 300)
sns.barplot(data=hippo_counts, x='cell_type', y='counts', hue='condition', palette='rocket')
plt.title('Hippocampus')
plt.xlabel('')
plt.ylabel('Number of Cells')
plt.xticks(rotation=45)  
plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(False) 
plt.show()

## Corpus Callosum

callosum_counts = callosum_data.groupby(['cell_type', 'condition'], observed=True).size().reset_index(name='counts')

plt.figure(figsize=(7, 5), dpi = 300)
sns.barplot(data=callosum_counts, x='cell_type', y='counts', hue='condition', palette='rocket')
plt.title('Corpus callosum')
plt.xlabel('')
plt.ylabel('Number of Cells')
plt.xticks(rotation=45)  
plt.legend(title='Group', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.grid(False) 
plt.show()

###############################################################################