#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 16:09:40 2024

@author: maureen
"""

import pandas as pd

def merge_annotations(master, subset, cell_type_col='Cell_type', subcluster_col='subcluster'):
    # Find common indices
    common_indices = master.obs.index.intersection(subset.obs.index)

    # 'Cell_type' in adata.obs
    if cell_type_col in subset.obs:
        # Ensure the master's column is categorical
        if not isinstance(master.obs[cell_type_col].dtype, pd.CategoricalDtype):
            master.obs[cell_type_col] = master.obs[cell_type_col].astype('category')
        # Ensure the subset's column is categorical
        if not isinstance(subset.obs[cell_type_col].dtype, pd.CategoricalDtype):
            subset.obs[cell_type_col] = subset.obs[cell_type_col].astype('category')

        master_categories = set(master.obs[cell_type_col].cat.categories)
        subset_categories = set(subset.obs[cell_type_col].cat.categories)
        new_categories = subset_categories - master_categories

        if new_categories:
            master.obs[cell_type_col] = master.obs[cell_type_col].cat.add_categories(new_categories)
        
        master.obs.loc[common_indices, cell_type_col] = subset.obs.loc[common_indices, cell_type_col].astype(master.obs[cell_type_col].dtype)

    # 'subcluster' in adata.obs
    if subcluster_col in subset.obs:
        # Ensure the master's column is categorical
        if not isinstance(master.obs[subcluster_col].dtype, pd.CategoricalDtype):
            master.obs[subcluster_col] = master.obs[subcluster_col].astype('category')
        # Ensure the subset's column is categorical
        if not isinstance(subset.obs[subcluster_col].dtype, pd.CategoricalDtype):
            subset.obs[subcluster_col] = subset.obs[subcluster_col].astype('category')

        master_categories = set(master.obs[subcluster_col].cat.categories)
        subset_categories = set(subset.obs[subcluster_col].cat.categories)
        new_categories = subset_categories - master_categories

        if new_categories:
            master.obs[subcluster_col] = master.obs[subcluster_col].cat.add_categories(new_categories)
        
        master.obs.loc[common_indices, subcluster_col] = subset.obs.loc[common_indices, subcluster_col].astype(master.obs[subcluster_col].dtype)

merge_annotations(adata, subset)