#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  5 15:10:09 2024

@author: maureen
"""

# Size factors function
def calculate_size_factors(subset_data):
    total_counts = subset_data.X.sum(axis=1)
    size_factors = total_counts / np.mean(total_counts)
    return size_factors

# DE function
def run_de_analysis(adata, cell_type, reference_group, comparison_group, min_cells_per_group=40, max_imbalance_ratio=0.9):
    subset_data = adata[(adata.obs['cell_type'] == cell_type) & 
                        (adata.obs['group'].isin([reference_group, comparison_group]))].copy()

    group_counts = subset_data.obs['group'].value_counts()
    if any(group_counts < min_cells_per_group):
        return f"Skipped: {reference_group} vs {comparison_group} (n is too low)"

    if max(group_counts) / sum(group_counts) > max_imbalance_ratio:
        return f"Skipped: {reference_group} vs {comparison_group} (imbalance)"

    subset_data.X = subset_data.layers['counts'].copy()

    sc.pp.filter_genes(subset_data, min_cells=3)

    if issparse(subset_data.X):
        subset_data.X = subset_data.X.toarray()

    size_factors = calculate_size_factors(subset_data)
    
    subset_data.obs['size_factors'] = size_factors

    test_result = de.test.wald(
        data=subset_data,
        formula_loc="~ 1 + group",
        factor_loc_totest="group",
        size_factors=size_factors
    )
    
    df_result = test_result.summary()
    df_result['cell_type'] = cell_type
    df_result['comparison'] = f"{reference_group}{comparison_group}"
    
    return df_result

###############################################################################

# Run differential expression

cell_types = ['Excitatory neuron', 'Inhibitory neuron', 'Mixed neuron', 'Astrocyte', 'Oligodendrocyte', 'OPC', 'Microglia', 'Choroid-plexus epithelial', 'Fibroblast']

comparisons = [('A', 'C'), ('C', 'E'), ('A', 'B'), ('A', 'D'), ('C', 'D'), ('D', 'F'), ('E', 'F')]

all_results = []

for cell_type in cell_types:
    for reference_group, comparison_group in comparisons:
        result = run_de_analysis(
            adata=adata,
            cell_type=cell_type,
            reference_group=reference_group,
            comparison_group=comparison_group
        )
        if isinstance(result, pd.DataFrame):
            all_results.append(result)
        else:
            print(f"Skipping {cell_type} {reference_group} vs {comparison_group}: {result}")  

final_results = pd.concat(all_results, ignore_index=True)

###############################################################################