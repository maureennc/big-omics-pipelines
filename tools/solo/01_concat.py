#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 08:10:22 2024

@author: maureen
"""

import os
import scanpy as sc
from scipy.sparse import csr_matrix
import anndata as ad

###############################################################################

# IMPORT DATA

data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/azenta-files/01_analysis/cellranger_count"

sample_A = sc.read_10x_mtx(os.path.join(data_dir, "1/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_B = sc.read_10x_mtx(os.path.join(data_dir, "3/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_C = sc.read_10x_mtx(os.path.join(data_dir, "5/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_D = sc.read_10x_mtx(os.path.join(data_dir, "6/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_E = sc.read_10x_mtx(os.path.join(data_dir, "7/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)
sample_F = sc.read_10x_mtx(os.path.join(data_dir, "8/filtered_feature_bc_matrix"), var_names='gene_symbols', make_unique=True)


## Create sample_ID column with original folder name
sample_A.obs['sample_id'] = '1'
sample_B.obs['sample_id'] = '3'
sample_C.obs['sample_id'] = '5'
sample_D.obs['sample_id'] = '6'
sample_E.obs['sample_id'] = '7'
sample_F.obs['sample_id'] = '8'

###############################################################################

# CONCATENATION

## Perform concatenation

adata_list = [sample_A, sample_B, sample_C, sample_D, sample_E, sample_F]

adata_concat = ad.concat(adata_list, 
                         label='group', 
                         keys=['A', 'B', 'C', 'D', 'E', 'F'], 
                         index_unique='-',
                         join = 'inner',
                         merge = 'unique') # prevents adata.var from being dropped

###############################################################################

# EXPORT 

save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/2/h5ad"

adata_concat.X = csr_matrix(adata_concat.X)
adata_concat.write_h5ad(os.path.join(save_dir, '1-raw-concat.h5ad'))

###############################################################################