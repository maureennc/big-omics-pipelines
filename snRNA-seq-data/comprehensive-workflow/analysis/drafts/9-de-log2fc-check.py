#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 17 10:23:56 2024

@author: maureen

- for all of them the values are bigger in group 2 (second alphabetically). So diffxpy sets group 1 alphabetically as reference
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc

################################################################################################################################

# IMPORT

## AnnData
data_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/h5ad"

adata = sc.read_h5ad(os.path.join(data_dir, '3-tbi-annotated-hvg.h5ad'))
adata.X = adata.layers['log1p'].copy()


################################################################################################################################

# READ CELL TYPE CSV FILES

csv_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de'

ac_results = pd.read_csv(os.path.join(csv_dir, 'cell_type/AC-comparison-cell_type-filtered.csv'))
ce_results = pd.read_csv(os.path.join(csv_dir, 'cell_type/CE-comparison-cell_type-filtered.csv'))
ab_results = pd.read_csv(os.path.join(csv_dir, 'cell_type/AB-comparison-cell_type-filtered.csv'))
ad_results = pd.read_csv(os.path.join(csv_dir, 'cell_type/AD-comparison-cell_type-filtered.csv'))
cd_results = pd.read_csv(os.path.join(csv_dir, 'cell_type/CD-comparison-cell_type-filtered.csv'))
df_results = pd.read_csv(os.path.join(csv_dir, 'cell_type/DF-comparison-cell_type-filtered.csv'))


################################################################################################################################

# CELL TYPE FUNCTION

def verify_foldchange_direction_cell_type(filtered_results, group1, group2, cell_types_list, adata):
    results = []

    for cell_type in cell_types_list:
        ## Filter the DataFrame for each cell type
        de_df = filtered_results[filtered_results['cell_type'] == cell_type]

        if de_df.empty:
            print(f"No results found for cell type: {cell_type} in comparison {group1} vs {group2}")
            continue

        ## Sort by log2fc to find the most upregulated genes, descending
        de_df_sorted = de_df.sort_values('log2fc', ascending=False).reset_index(drop=True)

        if de_df_sorted.empty:
            print(f"No sorted results found for cell type: {cell_type} in comparison {group1} vs {group2}")
            continue

        ## Get the most upregulated gene
        most_upregulated = de_df_sorted.iloc[0]['gene']
        log2fc = de_df_sorted.iloc[0]['log2fc']

        ## Get the index of the most upregulated gene in the AnnData var_names
        gene_index = np.where(adata.var_names == most_upregulated)[0]
        if gene_index.size == 0:
            print(f"Gene {most_upregulated} not found in AnnData for cell type: {cell_type}")
            continue
        i = gene_index[0]

        ## Extract expression data for the most upregulated gene from the two specified groups within the cell type
        filtered_adata = adata[(adata.obs['group'].isin([group1, group2])) & (adata.obs['cell_type'] == cell_type)]
        a = filtered_adata[filtered_adata.obs['group'] == group1, i].X
        b = filtered_adata[filtered_adata.obs['group'] == group2, i].X
        a_mean = np.mean(a)
        b_mean = np.mean(b)

        results.append({
            'cell_type': cell_type,
            'gene': most_upregulated,
            f'group_{group1}_mean': a_mean,
            f'group_{group2}_mean': b_mean,
            'log2FC': log2fc
        })

    return pd.DataFrame(results)

################################################################################################################################

# RUN AND EXPORT RESULTS

## Run summary for each comparison
cell_types_list = ['Astrocyte', 'Excitatory neuron', 'Inhibitory neuron', 'Oligodendrocyte', 'OPC', 'Microglia', 'Fibroblast', 'Unassigned']

ac_fc_check_1 = verify_foldchange_direction_cell_type(ac_results, 'A', 'C', cell_types_list, adata)
ce_fc_check_1 = verify_foldchange_direction_cell_type(ce_results, 'C', 'E', cell_types_list, adata)
ab_fc_check_1 = verify_foldchange_direction_cell_type(ab_results, 'A', 'B', cell_types_list, adata)
ad_fc_check_1 = verify_foldchange_direction_cell_type(ad_results, 'A', 'D', cell_types_list, adata)
cd_fc_check_1 = verify_foldchange_direction_cell_type(cd_results, 'C', 'D', cell_types_list, adata)
df_fc_check_1 = verify_foldchange_direction_cell_type(df_results, 'D', 'F', cell_types_list, adata)

#No results found for cell type: Fibroblast in comparison A vs D
#No results found for cell type: Microglia in comparison C vs D
#No results found for cell type: Fibroblast in comparison C vs D
#No results found for cell type: Fibroblast in comparison D vs F
#No results found for cell type: Unassigned in comparison D vs F

## Export
save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de"

ac_fc_check_1.to_csv(os.path.join(save_dir, 'cell_type/check/AC-top-upreg-cell_type.csv'), index = False)
ce_fc_check_1.to_csv(os.path.join(save_dir, 'cell_type/check/CE-top-upreg-cell_type.csv'), index = False)
ab_fc_check_1.to_csv(os.path.join(save_dir, 'cell_type/check/AB-top-upreg-cell_type.csv'), index = False)
ad_fc_check_1.to_csv(os.path.join(save_dir, 'cell_type/check/AD-top-upreg-cell_type.csv'), index = False)
cd_fc_check_1.to_csv(os.path.join(save_dir, 'cell_type/check/CD-top-upreg-cell_type.csv'), index = False)
df_fc_check_1.to_csv(os.path.join(save_dir, 'cell_type/check/DF-top-upreg-cell_type.csv'), index = False)

################################################################################################################################

# READ CLUSTER CSV FILES

csv_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de'

ac_results = pd.read_csv(os.path.join(csv_dir, 'cluster/AC-comparison-cluster-filtered.csv'))
ce_results = pd.read_csv(os.path.join(csv_dir, 'cluster/CE-comparison-cluster-filtered.csv'))
ab_results = pd.read_csv(os.path.join(csv_dir, 'cluster/AB-comparison-cluster-filtered.csv'))
ad_results = pd.read_csv(os.path.join(csv_dir, 'cluster/AD-comparison-cluster-filtered.csv'))
cd_results = pd.read_csv(os.path.join(csv_dir, 'cluster/CD-comparison-cluster-filtered.csv'))
df_results = pd.read_csv(os.path.join(csv_dir, 'cluster/DF-comparison-cluster-filtered.csv'))


################################################################################################################################

# CLUSTER FUNCTION

def verify_foldchange_direction_cluster(filtered_results, group1, group2, clusters_list, adata):
    results = []

    for cluster in clusters_list:
        ## Filter the DataFrame for each cell type
        de_df = filtered_results[filtered_results['cluster'] == cluster]

        if de_df.empty:
            print(f"No results found for cluster: {cluster} in comparison {group1} vs {group2}")
            continue

        ## Sort by log2fc to find the most upregulated genes, descending
        de_df_sorted = de_df.sort_values('log2fc', ascending=False).reset_index(drop=True)

        if de_df_sorted.empty:
            print(f"No sorted results found for cluster: {cluster} in comparison {group1} vs {group2}")
            continue

        ## Get the most upregulated gene
        most_upregulated = de_df_sorted.iloc[0]['gene']
        log2fc = de_df_sorted.iloc[0]['log2fc']

        ## Get the index of the most upregulated gene in the AnnData var_names
        gene_index = np.where(adata.var_names == most_upregulated)[0]
        if gene_index.size == 0:
            print(f"Gene {most_upregulated} not found in AnnData for cell type: {cluster}")
            continue
        i = gene_index[0]

        ## Extract expression data for the most upregulated gene from the two specified groups within the cell type
        filtered_adata = adata[(adata.obs['group'].isin([group1, group2])) & (adata.obs['cluster'] == cluster)]
        a = filtered_adata[filtered_adata.obs['group'] == group1, i].X
        b = filtered_adata[filtered_adata.obs['group'] == group2, i].X
        a_mean = np.mean(a)
        b_mean = np.mean(b)

        results.append({
            'cluster': cluster,
            'gene': most_upregulated,
            f'group_{group1}_mean': a_mean,
            f'group_{group2}_mean': b_mean,
            'log2FC': log2fc
        })

    return pd.DataFrame(results)

################################################################################################################################

# RUN AND EXPORT RESULTS

## Run summary for each comparison
clusters_list = adata.obs['cluster'].unique().tolist()

ac_fc_check_2 = verify_foldchange_direction_cluster(ac_results, 'A', 'C', clusters_list, adata)
ce_fc_check_2 = verify_foldchange_direction_cluster(ce_results, 'C', 'E', clusters_list, adata)
ab_fc_check_2 = verify_foldchange_direction_cluster(ab_results, 'A', 'B', clusters_list, adata)
ad_fc_check_2 = verify_foldchange_direction_cluster(ad_results, 'A', 'D', clusters_list, adata)
cd_fc_check_2 = verify_foldchange_direction_cluster(cd_results, 'C', 'D', clusters_list, adata)
df_fc_check_2 = verify_foldchange_direction_cluster(df_results, 'D', 'F', clusters_list, adata)


## Export
save_dir = "/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/results/spreadsheets/de"

ac_fc_check_2.to_csv(os.path.join(save_dir, 'cluster/check/AC-top-upreg-cluster.csv'), index = False)
ce_fc_check_2.to_csv(os.path.join(save_dir, 'cluster/check/CE-top-upreg-cluster.csv'), index = False)
ab_fc_check_2.to_csv(os.path.join(save_dir, 'cluster/check/AB-top-upreg-cluster.csv'), index = False)
ad_fc_check_2.to_csv(os.path.join(save_dir, 'cluster/check/AD-top-upreg-cluster.csv'), index = False)
cd_fc_check_2.to_csv(os.path.join(save_dir, 'cluster/check/CD-top-upreg-cluster.csv'), index = False)
df_fc_check_2.to_csv(os.path.join(save_dir, 'cluster/check/DF-top-upreg-cluster.csv'), index = False)


################################################################################################################################
