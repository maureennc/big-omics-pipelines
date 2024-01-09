#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan  9 13:54:08 2024

@author: maureen
"""

# SCRUBLET SINGLE SAMPLE

def run_scrublet_single(adata):
    # Run Scrublet
    scrub = scr.Scrublet(adata.X, expected_doublet_rate=0.20)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                              min_cells=3, 
                                                              min_gene_variability_pctl=85, 
                                                              n_prin_comps=30)
    adata.obs['doublet_scores'] = doublet_scores
    adata.obs['predicted_doublets'] = predicted_doublets

    # Plot Histograms
    plt.figure(figsize=(10, 6))
    sns.histplot(scrub.doublet_scores_obs_, bins=30, color="blue", label="Observed", kde=True)
    sns.histplot(scrub.doublet_scores_sim_, bins=30, color="red", label="Simulated", kde=True)
    plt.xlabel('Doublet Score')
    plt.ylabel('Density')
    plt.legend()
    plt.show()

    # Extract cell barcodes and store Scrublet data
    scrublet_rows = [{
        'Cell_Barcode': barcode,
        'Observed_Score': obs_score, 
        'Simulated_Score': sim_score, 
        'Predicted_Doublet': pred_doublet
    } for barcode, obs_score, sim_score, pred_doublet in zip(adata.obs.index, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_, predicted_doublets)]

    # Create a DataFrame from the list of rows
    scrublet_df = pd.DataFrame(scrublet_rows)
    return scrublet_df
