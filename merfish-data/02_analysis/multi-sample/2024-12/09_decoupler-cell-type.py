#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 08:37:51 2024

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
import squidpy as sq
import random

###############################################################################

# SETTINGS


## Matplotlib
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['figure.dpi'] = 500
random.seed(0)

###############################################################################

# IMPORT DATA

data_dir = '/Users/maureen/Desktop/toxo-merfish/data/h5ad/processed'

adata = sc.read_h5ad(os.path.join(data_dir, 'concat-split-trained-annotated.h5ad'))

###############################################################################

# SPATIAL PROCESSING

sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key='cell_type')
sq.pl.nhood_enrichment(adata, cluster_key='cell_type', figsize=(2, 2), palette = 'rocket', cmap = 'rocket')


sq.gr.spatial_neighbors(E009)
sq.gr.nhood_enrichment(E009, cluster_key='recluster')
sq.pl.nhood_enrichment(E009, cluster_key='recluster', figsize=(2, 2))

