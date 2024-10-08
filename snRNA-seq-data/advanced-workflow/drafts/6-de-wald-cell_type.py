#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun 16 12:26:54 2024

@author: maureen

DE testing criteria
- Total cell type sample size > 100 cells 
- Individual group sample size > 40 cells
- Group balance; no group accounts for less than 10% of the total number of cells in the cell type

"""

import os
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse
import diffxpy.api as de

################################################################################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '3-tbi-annotated-hvg.h5ad'))

################################################################################################################################

# PREPARE ANNDATA

## Mask artifact clusters
adata = adata[adata.obs['cell_class'] != 'Artifact'].copy()

## Restore the raw counts for updated normalization
adata.X = adata.layers['log1p'].copy()

## Convert to dense array
if issparse(adata.X):
    adata.X = adata.X.toarray()

################################################################################################################################

# CREATE ANALYSIS SUBSETS

ac_data = adata[adata.obs['group'].isin(['A', 'C'])].copy()
ce_data = adata[adata.obs['group'].isin(['C', 'E'])].copy()
ab_data = adata[adata.obs['group'].isin(['A', 'B'])].copy()

ad_data = adata[adata.obs['group'].isin(['A', 'D'])].copy()
cd_data = adata[adata.obs['group'].isin(['C', 'D'])].copy()
df_data = adata[adata.obs['group'].isin(['D', 'F'])].copy()

################################################################################################################################

# AC - DE

ac_results_list = []
cell_types = ac_data.obs['cell_type'].unique()

## Loop through each  cell type and prepare the data
for cell_type in cell_types:
    cell_type_data = ac_data[ac_data.obs['cell_type'] == cell_type].copy()
    groups = cell_type_data.obs['group'].value_counts()
    
    ## Filter out groups with less than 40 cells for specific cell type comparison
    groups = groups[groups >= 40].index.tolist()
    
    ## Requirements for total cell type sample size > 100 cells and balance (> 10% representation)
    if cell_type_data.shape[0] > 100 and len(groups) > 1:
        balanced = all(cell_type_data.obs['group'].value_counts(normalize=True) > 0.10)
        
        if balanced:
            try:
                ## Perform DE
                test_result = de.test.wald(
                    data=cell_type_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cell_type'] = cell_type
                ac_results_list.append(df_result)
                print(f"Processed cell type: {cell_type}")
            except Exception as e:
                print(f"Error processing cell type {cell_type}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cell_type} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cell_type} due to insufficient data or group criteria not met.")

## Concatenate
if ac_results_list:
    results_cell_type = pd.concat(ac_results_list, ignore_index=True)
    print(results_cell_type)
else:
    print("No results to display.")
    
ac_results = pd.concat(ac_results_list, ignore_index=True)

################################################################################################################################

# CE - CELL TYPE DE

ce_results_list = []
cell_types = ce_data.obs['cell_type'].unique()

## Loop through each  cell type and prepare the data
for cell_type in cell_types:
    cell_type_data = ce_data[ce_data.obs['cell_type'] == cell_type].copy()
    groups = cell_type_data.obs['group'].value_counts()
    
    ## Filter out groups with less than 40 cells for specific cell type comparison
    groups = groups[groups >= 40].index.tolist()
    
    ## Requirements for total cell type sample size > 100 cells and balance (> 10% representation)
    if cell_type_data.shape[0] > 100 and len(groups) > 1:
        balanced = all(cell_type_data.obs['group'].value_counts(normalize=True) > 0.10)
        
        if balanced:
            try:
                ## Perform DE
                test_result = de.test.wald(
                    data=cell_type_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cell_type'] = cell_type
                ce_results_list.append(df_result)
                print(f"Processed cell type: {cell_type}")
            except Exception as e:
                print(f"Error processing cell type {cell_type}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cell_type} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cell_type} due to insufficient data or group criteria not met.")

## Concatenate
if ce_results_list:
    results_cell_type = pd.concat(ce_results_list, ignore_index=True)
    print(results_cell_type)
else:
    print("No results to display.")
    
ce_results = pd.concat(ce_results_list, ignore_index=True)

################################################################################################################################

# AB - CELL TYPE DE

ab_results_list = []
cell_types = ab_data.obs['cell_type'].unique()

## Loop through each  cell type and prepare the data
for cell_type in cell_types:
    cell_type_data = ab_data[ab_data.obs['cell_type'] == cell_type].copy()
    groups = cell_type_data.obs['group'].value_counts()
    
    ## Filter out groups with less than 40 cells for specific cell type comparison
    groups = groups[groups >= 40].index.tolist()
    
    ## Requirements for total cell type sample size > 100 cells and balance (> 10% representation)
    if cell_type_data.shape[0] > 100 and len(groups) > 1:
        balanced = all(cell_type_data.obs['group'].value_counts(normalize=True) > 0.10)
        
        if balanced:
            try:
                ## Perform DE
                test_result = de.test.wald(
                    data=cell_type_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cell_type'] = cell_type
                ab_results_list.append(df_result)
                print(f"Processed cell type: {cell_type}")
            except Exception as e:
                print(f"Error processing cell type {cell_type}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cell_type} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cell_type} due to insufficient data or group criteria not met.")

## Concatenate
if ab_results_list:
    results_cell_type = pd.concat(ab_results_list, ignore_index=True)
    print(results_cell_type)
else:
    print("No results to display.")
    
ab_results = pd.concat(ab_results_list, ignore_index=True)


################################################################################################################################

# AD COMPARISON

ad_results_list = []
cell_types = ad_data.obs['cell_type'].unique()

## Loop through each  cell type and prepare the data
for cell_type in cell_types:
    cell_type_data = ad_data[ad_data.obs['cell_type'] == cell_type].copy()
    groups = cell_type_data.obs['group'].value_counts()
    
    ## Filter out groups with less than 40 cells for specific cell type comparison
    groups = groups[groups >= 40].index.tolist()
    
    ## Requirements for total cell type sample size > 100 cells and balance (> 10% representation)
    if cell_type_data.shape[0] > 100 and len(groups) > 1:
        balanced = all(cell_type_data.obs['group'].value_counts(normalize=True) > 0.10)
        
        if balanced:
            try:
                ## Perform DE
                test_result = de.test.wald(
                    data=cell_type_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cell_type'] = cell_type
                ad_results_list.append(df_result)
                print(f"Processed cell type: {cell_type}")
            except Exception as e:
                print(f"Error processing cell type {cell_type}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cell_type} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cell_type} due to insufficient data or group criteria not met.")

## Concatenate
if ad_results_list:
    results_cell_type = pd.concat(ad_results_list, ignore_index=True)
    print(results_cell_type)
else:
    print("No results to display.")
    
ad_results = pd.concat(ad_results_list, ignore_index=True)

################################################################################################################################

# CD

cd_results_list = []
cell_types = cd_data.obs['cell_type'].unique()

## Loop through each  cell type and prepare the data
for cell_type in cell_types:
    cell_type_data = cd_data[cd_data.obs['cell_type'] == cell_type].copy()
    groups = cell_type_data.obs['group'].value_counts()
    
    ## Filter out groups with less than 40 cells for specific cell type comparison
    groups = groups[groups >= 40].index.tolist()
    
    ## Requirements for total cell type sample size > 100 cells and balance (> 10% representation)
    if cell_type_data.shape[0] > 100 and len(groups) > 1:
        ## Check for balance: no group should be less than 10% of the cell type size
        balanced = all(cell_type_data.obs['group'].value_counts(normalize=True) > 0.10)
        
        if balanced:
            try:
                ## Perform DE
                test_result = de.test.wald(
                    data=cell_type_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cell_type'] = cell_type
                cd_results_list.append(df_result)
                print(f"Processed cell type: {cell_type}")
            except Exception as e:
                print(f"Error processing cell type {cell_type}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cell_type} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cell_type} due to insufficient data or group criteria not met.")

## Concatenate
if cd_results_list:
    results_cell_type = pd.concat(cd_results_list, ignore_index=True)
    print(results_cell_type)
else:
    print("No results to display.")
    
cd_results = pd.concat(cd_results_list, ignore_index=True)

################################################################################################################################

# DF COMPARISON

df_results_list = []
cell_types = df_data.obs['cell_type'].unique()

## Loop through each  cell type and prepare the data
for cell_type in cell_types:
    cell_type_data = df_data[df_data.obs['cell_type'] == cell_type].copy()
    groups = cell_type_data.obs['group'].value_counts()
    
    ## Filter out groups with less than 40 cells for specific cell type comparison
    groups = groups[groups >= 40].index.tolist()
    
    ## Requirements for total cell type sample size > 100 cells and balance (> 10% representation)
    if cell_type_data.shape[0] > 100 and len(groups) > 1:
        ## Check for balance: no group should be less than 10% of the cell type size
        balanced = all(cell_type_data.obs['group'].value_counts(normalize=True) > 0.10)
        
        if balanced:
            try:
                ## Perform DE
                test_result = de.test.wald(
                    data=cell_type_data,
                    formula_loc="~ 1 + group",
                    factor_loc_totest="group"
                )
                df_result = test_result.summary()
                df_result['cell_type'] = cell_type
                df_results_list.append(df_result)
                print(f"Processed cell type: {cell_type}")
            except Exception as e:
                print(f"Error processing cell type {cell_type}: {e}")
        else:
            print(f"Skipping DE analysis for cell type {cell_type} due to imbalance in group sizes.")
    else:
        print(f"Skipping DE analysis for cell type {cell_type} due to insufficient data or group criteria not met.")

## Concatenate
if df_results_list:
    results_cell_type = pd.concat(df_results_list, ignore_index=True)
    print(results_cell_type)
else:
    print("No results to display.")
    
df_results = pd.concat(df_results_list, ignore_index=True)

################################################################################################################################

# FILTER EXTREME LOG2FC WITH WEAK EXPRESSION
## This occurs with diffxpy when there are no reads for a gene detected in one group due to very low expression
ac_filtered = ac_results[(abs(ac_results['log2fc']) <= 100) & (ac_results['mean'] >= 0.05)].copy()
ce_filtered = ce_results[(abs(ce_results['log2fc']) <= 100) & (ce_results['mean'] >= 0.05)].copy()
ab_filtered = ab_results[(abs(ab_results['log2fc']) <= 100) & (ab_results['mean'] >= 0.05)].copy()

ad_filtered = ad_results[(abs(ad_results['log2fc']) <= 100) & (ad_results['mean'] >= 0.05)].copy()
cd_filtered = cd_results[(abs(cd_results['log2fc']) <= 100) & (cd_results['mean'] >= 0.05)].copy()
df_filtered = df_results[(abs(df_results['log2fc']) <= 100) & (df_results['mean'] >= 0.05)].copy()

################################################################################################################################

# EXPORT CSV FILES

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de/cell_type'

ac_results.to_csv(os.path.join(save_dir, 'full/AC-comparison-cell_type-full.csv'), index = False)
ce_results.to_csv(os.path.join(save_dir, 'full/CE-comparison-cell_type-full.csv'), index = False)
ab_results.to_csv(os.path.join(save_dir, 'full/AB-comparison-cell_type-full.csv'), index = False)
ad_results.to_csv(os.path.join(save_dir, 'full/AD-comparison-cell_type-full.csv'), index = False)
cd_results.to_csv(os.path.join(save_dir, 'full/CD-comparison-cell_type-full.csv'), index = False)
df_results.to_csv(os.path.join(save_dir, 'full/DF-comparison-cell_type-full.csv'), index = False)

ac_filtered.to_csv(os.path.join(save_dir, 'filtered/AC-comparison-cell_type-filtered.csv'), index = False)
ce_filtered.to_csv(os.path.join(save_dir, 'filtered/CE-comparison-cell_type-filtered.csv'), index = False)
ab_filtered.to_csv(os.path.join(save_dir, 'filtered/AB-comparison-cell_type-filtered.csv'), index = False)
ad_filtered.to_csv(os.path.join(save_dir, 'filtered/AD-comparison-cell_type-filtered.csv'), index = False)
cd_filtered.to_csv(os.path.join(save_dir, 'filtered/CD-comparison-cell_type-filtered.csv'), index = False)
df_filtered.to_csv(os.path.join(save_dir, 'filtered/DF-comparison-cell_type-filtered.csv'), index = False)

################################################################################################################################

#filtered_de_cell_type = results_cell_type[(abs(results_cell_type['log2fc']) <= 100) & (results_cell_type['mean'] >= 0.05)]
#filtered_de_cell_type['log2fc'] = -filtered_de_cell_type['log2fc']
#filtered_de_cell_type.to_csv(os.path.join(save_dir, 'mg-sting-wt-ko-DE-cell_type-median-norm.csv'), index = False)
