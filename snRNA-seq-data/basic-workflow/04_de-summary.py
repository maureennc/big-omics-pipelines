#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 25 17:00:04 2024

@author: maureen
"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats
import random


################################################################################################################################

# SETTINGS

## Random seed
random.seed(0)
np.random.seed(0)

## Matplotlib
%matplotlib qt5
plt.rcParams['font.family'] = 'Arial'


################################################################################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/mg-sting-ko/results/spreadsheets/de"

filtered_de_cell_type = pd.read_csv(os.path.join(data_dir, 'mg-sting-wt-ko-DE-cell_type.csv'))
filtered_de_cluster = pd.read_csv(os.path.join(data_dir, 'mg-sting-wt-ko-DE-cluster.csv'))

################################################################################################################################

de_criteria = (filtered_de_cluster['qval'] < 0.05) & (abs(filtered_de_cluster['log2fc']) >= 0.5)

# PLOT CELL TYPE
plt.figure(figsize=(3, 3), dpi=300)
de_cell_type_counts = filtered_de_cell_type[de_criteria].groupby('cell_type').size().reset_index(name='counts')
sns.barplot(data=de_cell_type_counts, x='cell_type', y='counts', palette = 'viridis')
plt.title('Number of DEGs by Cell Type')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


# PLOT CLUSTER
plt.figure(figsize=(6, 3), dpi=300)
de_cluster_counts = filtered_de_cluster[de_criteria].groupby('cluster').size().reset_index(name='counts')
sns.barplot(data=de_cluster_counts, x='cluster', y='counts', palette = 'viridis')
plt.title('Number of DEGs by Cluster')
plt.xticks(rotation=45, ha='right')
plt.tight_layout()
plt.show()


################################################################################################################################

## Save directory
#save_dir = "/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-3/nestin-sting-ko/results/spreadsheets/de"

