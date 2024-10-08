#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 16 11:34:49 2023

@author: maureen
"""

import os
os.chdir("/Users/maureen/Documents/projects/harris-lab/Lydia/MERFISH_data/2023_09_08_Brain-28DPI-first_fresh_frozen")

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc
sc.set_figure_params(scanpy = True, dpi = 150, dpi_save = 400)

import squidpy as sq
import scrublet as scr

# environment = 'sc'

################################

# Import data and construct adata objects

vizgen_dir_diameter = Path().resolve() / "vpt_cyto2A-C_diameter" / "squidpy"


adata = sq.read.vizgen(path = vizgen_dir_diameter, 
                       counts_file = "2B_cell_by_gene.csv",
                       meta_file = "2B_cell_metadata.csv",
                       transformation_file = "micron_to_mosaic_pixel_transform.csv")

#################################

anndata_list = [adata]

def pp_and_filter(anndata_list):
    for adata in anndata_list:
        # Standard filter for genes and cells
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.filter_cells(adata, min_genes=10)
        sc.pp.filter_cells(adata, min_counts=50)
        sc.pp.filter_cells(adata, max_counts=3000)
        
        # Volume filter
        min_volume = 50
        max_volume = 4000

        adata._inplace_subset_obs(adata.obs['volume'] > min_volume)
        adata._inplace_subset_obs(adata.obs['volume'] < max_volume)
        
        # QC metric calculation
        sc.pp.calculate_qc_metrics(adata, percent_top=(50, 100, 200, 300), inplace=True)

# Process and filter each Anndata object in-place
pp_and_filter(anndata_list)



###############

# Scrublet

def scrub_init(anndata_list, expected_doublet_rate=0.5):
    for i, adata in enumerate(anndata_list):
        # Scrublet step
        scrub = scr.Scrublet(adata.X, expected_doublet_rate=expected_doublet_rate)
        adata.obs['doublet_scores'], adata.obs['predicted_doublets'] = scrub.scrub_doublets(
            min_counts=10,
            min_cells=3,
            min_gene_variability_pctl=85
        )
        
        # Scrublet histogram plot for the current dataset
        plt.figure(figsize=(8, 4))
        scrub.plot_histogram()
        plt.savefig(f'scrublet_histogram_{i}.png')
        plt.close()

scrub_init(anndata_list)

#################
# create backup objects before proceeding

adata_backup = adata

###########################

# Manually exclude doublets based on scrublet score
adata = adata[adata.obs["doublet_scores"] < 0.60]


######################################

# CLUSTERING
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=20)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5)
sc.pl.umap(adata, color = ['leiden'])


##########################################

# CLUSTER ANNOTATION

sc.tl.rank_genes_groups(adata, 'leiden')
sc.pl.rank_genes_groups(adata, n_genes=20, sharey=False)

Cell_type = {
"0": "Excitatory neuron 1",
"1": "Astrocyte",
"2": "Inhibitory neuron",
"3": "Microglia",
"4": "Oligodendrocyte",
"5": "Endothelial",
"6": "Macrophages",
"7": "Excitatory neuron 2",
"8": "CD8+ T cell",
"9": "Excitatory neuron 3",
"10": "CD4+ T cell",
"11": "Pericyte",
"12": "Excitatory neuron 4",
"13": "OPC",
"14": "Excitatory neuron 5"
}


category= {
"0": "Excitatory neuron",
"1": "Astrocyte",
"2": "Inhibitory neuron",
"3": "Microglia",
"4": "Oligodendrocyte",
"5": "Endothelial",
"6": "Macrophages",
"7": "Excitatory neuron",
"8": "CD8+ T cell",
"9": "Excitatory neuron",
"10": "CD4+ T cell",
"11": "Pericyte",
"12": "Excitatory neuron",
"13": "OPC",
"14": "Excitatory neuron"
}

adata.obs['Cell_type'] = adata.obs.leiden.map(Cell_type)
adata.obs['category'] = adata.obs.leiden.map(category)


sc.pl.umap(adata, color = ['Cell_type'])

##########################################

# COMPOSITIONAL ANALYSIS
cell_counts = adata.obs['Cell_type'].value_counts()
cell_counts = pd.DataFrame(cell_counts)
cell_counts['Cell_type'] = cell_counts.index
sns.barplot(cell_counts, x = 'count', y = 'Cell_type')


cell_counts2 = adata.obs['category'].value_counts()
cell_counts2 = pd.DataFrame(cell_counts2)
cell_counts2['category'] = cell_counts2.index
sns.barplot(cell_counts2, x = 'count', y = 'category')


##########################################

# HEATMAP

sc.tl.rank_genes_groups(adata, 'category')
sc.pl.rank_genes_groups_heatmap(adata, n_genes = 5, show_gene_labels=True)


# DOTPLOT

gene_list = ['Rbfox3', 'Cx3cl1', 'Slc17a7', 'Gad2', 'Atp1a2', 'Aqp4', 'Aldh1l1', 'Sox10', 'Olig1', 'Gjb1', 'Pdgfra', 'Cd3e', 'Cd4', 'Cd8a', 'Mpeg1', 'Csf1r', 'Nos2', 'Tmem119',  'Cldn5', 'Pecam1', 'Cdh5', 'Rgs5', 'Acta2']
sc.pl.dotplot(adata, gene_list, groupby = ['category'], dendrogram = True)

##########################################


# SPATIAL SCATTERPLOT BY CELL TYPE

sq.gr.spatial_neighbors(adata, key_added = 'spatial')

sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Excitatory neuron'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Inhibitory neuron'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Oligodendrocyte'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['OPC'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Astrocyte'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Endothelial'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Pericyte'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['CD4+ T cell'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['CD8+ T cell'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Microglia'], wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['Macrophages'], wspace=0.4,)

sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['CD8+ T cell', 'CD4+ T cell'], wspace=0.4,)

##########################################


# MORAN'S I
sc.pp.highly_variable_genes(adata, n_top_genes = 338)

genes = adata[:, adata.var.highly_variable].var_names.values[:338]

sq.gr.spatial_neighbors(adata)
sq.gr.spatial_autocorr(
    adata,
    mode="moran",
    genes=genes,
    n_perms=100,
    n_jobs=1,
)
adata.uns["moranI"].head(30)
adata.uns["moranI"].tail(30)

moranI = pd.DataFrame(adata.uns['moranI'])

sq.pl.spatial_scatter(adata, shape=None,color=['Slc17a7'],cmap='magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Gas6'],cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Rora'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Ifng'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Tnf'], cmap = 'magma', wspace=0.4,)

sq.pl.spatial_scatter(adata, shape=None,color=['C3'],cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Rbfox3'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Cxcl10'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Sox10'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Cx3cl1'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Tap1'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Il33'], cmap = 'magma', wspace=0.4,)
sq.pl.spatial_scatter(adata, shape=None,color=['Stat1'], cmap = 'magma', wspace=0.4,)



##########################################

# GENE SCORES

## print modules:

#import csv

#def read_genes_from_csv(file_path):
#    gene_list = []
#    with open(file_path, 'r') as csvfile:
#        reader = csv.reader(csvfile)
#        for row in reader:
#            gene_list.append(row[0])
#    return gene_list

#def process_csv_files(folder_path):
 #   gene_lists = []
  #  for filename in os.listdir(folder_path):
   #     if filename.endswith(".csv"):
    #        file_path = os.path.join(folder_path, filename)
     #       genes = read_genes_from_csv(file_path)
      #      gene_lists.append(genes)
    #return gene_lists

# Replace 'your_folder_path' with the actual path to the folder containing your CSV files
#folder_path = '/Users/maureen/Documents/projects/harris-lab/Tajie/TH_faculty_meeting_2023-12-1/gene_modules'
#all_gene_lists = process_csv_files(folder_path)

# Now 'all_gene_lists' is a list of lists, where each inner list contains gene names from a CSV file
#for i, genes in enumerate(all_gene_lists):
 #   print(f"Genes from file {i+1}: {genes}")


## create list
cytokine = ['Il1a', 'Il1b', 'Il1r1', 'Il1r2', 'Il2', 'Il2ra', 'Il2rb', 'Il2rg', 'Il4', 'Il4ra', 'Il6', 'Il6ra', 'Il7', 'Il7r', 'Il9', 'Il9r', 'Il10', 'Il10ra', 'Il10rb', 'Il12a', 'Il12b', 'Il12rb1', 'Il12rb2', 'Il13', 'Il13ra1', 'Il15', 'Il15ra', 'Il17a', 'Il17b', 'Il17c', 'Il17d', 'Il17f', 'Il17ra', 'Il17rb', 'Il17rc', 'Il18', 'Il21', 'Il21r', 'Il23a', 'Il23r', 'Il27', 'Il27ra', 'Il6st', 'Il33', 'Il1rl1', 'Il1rap', 'Il34', 'Csf1r', 'Ebi3', 'Il1rl2', 'Ifnab', 'Ifnb1', 'Ifnar1', 'Ifnar2', 'Ifng', 'Ifngr1', 'Ifngr2', 'Ifnlr1', 'Myd88', 'Socs1']
chemokine = ['Ccl2', 'Ccl19', 'Ccr7', 'Ccr5', 'Cxcl1', 'Cxcl2', 'Cxcr2', 'Cxcl9', 'Cxcl10', 'Cxcr3', 'Cx3cl1', 'Cxcr6', 'Cxcl16', 'Xcr1']
gf = ['Csf1', 'Csf2ra', 'Csf2rb', 'Csf2rb2', 'Csf3', 'Csf3r', 'Vegfa', 'Vegfb', 'Vegfc', 'Vegfd', 'Flt4', 'Tgfb1', 'Tgfb2', 'Tgfb3', 'Tgfbr1', 'Bdnf', 'Ntrk2', 'Ngf', 'Ntrk1']
cell_activation = ['Itgal', 'Cd44', 'Itga1', 'Itga4', 'B2m', 'Vcam1', 'Icam1', 'Prf1', 'Tap1', 'Mmp9', 'Cd69', 'Cd80', 'Cd86', 'Klrg1', 'Fos', 'Fosb', 'Jun', 'Nr4a1', 'Myc', 'Egr1', 'Arc', 'Nfatc1', 'Nfatc2', 'Nfatc3', 'Nfatc4', 'Nfat5', 'Nfkb1', 'Nfkb2', 'Rel', 'Rela', 'Relb', 'Stat1', 'Stat2', 'Stat3', 'Stat4', 'Stat5a', 'Stat5b', 'Stat6']
cell_death = ['Fasl', 'Ager', 'Casp8', 'Casp3', 'Casp7', 'Casp9', 'Tnfrsf1a', 'Tnfrsf1b', 'Ripk1', 'Ripk3', 'Fas', 'Tnfrsf10b', 'Tnfrsf25', 'Tnfrsf21', 'Fadd', 'Tradd', 'Apaf1', 'Cflar', 'Birc2', 'Birc3', 'Mlkl', 'Cyld', 'Traf2', 'Chuk', 'Ikbkg', 'Nfkbia', 'Tlr4', 'Ticam1', 'Zbp1', 'Bax', 'Bak1', 'Bad', 'Bid', 'Bcl2l11', 'Pmaip1', 'Bbc3', 'Bmf', 'Blk', 'Hrk', 'Casp1', 'Casp4', 'Nlrp3', 'Nlrc4', 'Nlrc3', 'Nlrp1a', 'Nlrp1b', 'Nlrp12', 'Naip6', 'Naip2', 'Nlrp1c-ps', 'Aim2', 'Gsdma', 'Gsdmd', 'Gsdme', 'Dpp8', 'Dpp9', 'Nek7', 'Pycard', 'S100a9', 'Tnf']
dam_activation = ['Axl', 'Msr1', 'Itgax', 'Spp1', 'Arg1', 'Siglec1', 'Gpnmb', 'Ccrl2', 'Lpl', 'Trem2', 'Clec7a', 'Syk', 'Mpeg1']

              
## score
sc.tl.score_genes(adata, gene_list = cytokine)
adata.obs.rename(columns={'score': 'cytokine_score'}, inplace=True)

sc.tl.score_genes(adata, chemokine)
adata.obs.rename(columns={'score': 'chemokine_score'}, inplace=True)

sc.tl.score_genes(adata, gf)
adata.obs.rename(columns={'score': 'gf_signaling'}, inplace=True)

sc.tl.score_genes(adata, cell_activation)
adata.obs.rename(columns={'score': 'cell_activation_score'}, inplace=True)

sc.tl.score_genes(adata, dam_activation)
adata.obs.rename(columns={'score': 'dam_score'}, inplace=True)

sc.tl.score_genes(adata, cell_death)
adata.obs.rename(columns={'score': 'cell_death_score'}, inplace=True)

## Visualize gene scores

sc.pl.matrixplot(adata, 
                 var_names = ['cytokine_score', 'chemokine_score', 'gf_signaling', 'cell_activation_score', 'dam_score', 'cell_death_score'], 
                 groupby='category', 
                 dendrogram=True,
                 standard_scale = 'var')


sc.pl.violin(adata, 
             keys='cytokine_score', 
             groupby='category', 
             rotation=90,
             aspect = 0.5)


##########################################


# OLIGO ACTIVATION

oligo_activation = ['Olig1', 'Gjb1', 'Stat1', 'Stat2', 'Stat3', 'Tap1', 'Il33', 'Il1rap', 'Il12rb1', 'Csf1', 'Zbp1', 'Jun', 'Arc']

## umaps
sc.pl.umap(adata, color = ['category'])
sc.pl.umap(adata, color = ['Stat1'])
sc.pl.umap(adata, color = ['Stat3'])
sc.pl.umap(adata, color = ['Tap1'])
sc.pl.umap(adata, color = ['Il1rap'])
sc.pl.umap(adata, color = ['Il12rb1'])
sc.pl.umap(adata, color = ['Arc'])
sc.pl.umap(adata, color = ['Il33'])


## violin plot
sc.pl.dotplot(adata, oligo_activation, groupby = 'category', dendrogram = True)


##########################################

# EXPORT

adata2 = adata
adata.uns['spatial'] #transformation matrix can't be a dictionary
del adata2.uns['spatial']

adata.write_h5ad('/Users/maureen/Documents/projects/harris-lab/Tajie/TH_faculty_meeting_2023-12-1/adata_clustered.hdf5')

adata = adata2
##########################################



# READ

adata = sc.read_h5ad('/Users/maureen/Documents/projects/harris-lab/Tajie/TH_faculty_meeting_2023-12-1/adata_clustered.hdf5')


##########################################


# GENERAL PLOTS

sc.pl.umap(adata, color = 'leiden')
sc.pl.umap(adata, color = 'category')

sc.pl.umap(adata, color = 'Rbfox3')
sc.pl.umap(adata, color = 'Nefm')
sc.pl.umap(adata, color = 'Ptprc')

sc.pl.umap(adata, color = 'Itgam')
sc.pl.umap(adata, color = 'Csf1r')
sc.pl.umap(adata, color = 'Cd3e')
sc.pl.umap(adata, color = 'Olig1')
sc.pl.umap(adata, color = 'Pecam1')



sq.pl.spatial_scatter(adata, shape=None,color=['category'], groups = ['leiden'], wspace=0.4,)


##########################################

# Immune recluster

immune = adata[adata.obs['leiden'].isin(['3', '6', '8', '10'])]

sc.tl.pca(immune, svd_solver='arpack')
sc.pp.neighbors(immune, n_neighbors=10)
sc.tl.umap(immune)
sc.pl.umap(immune, color = ['leiden'])

sc.tl.leiden(immune, resolution=.8, restrict_to=('leiden', ['3', '6', '8', '10']))


sc.pl.umap(immune, color = ['leiden_R', 'Csf1r', 'P2ry12', 'Mrc1'])
          

# rank genes
sc.tl.rank_genes_groups(immune, 'leiden_R')
sc.pl.rank_genes_groups(immune, n_genes=20, sharey=False)


# Rename clusters

immune_dict = {
    '3-6-8-10,0': 'Inf_mac_1',
    '3-6-8-10,1': 'Cd4',
    '3-6-8-10,2': 'Mg_1',
    '3-6-8-10,3': 'Cd8_1',
    '3-6-8-10,4': 'Mg_2',
    '3-6-8-10,5': 'Mg_3',
    '3-6-8-10,6': 'Mg_4',
    '3-6-8-10,7': 'Mg_5',
    '3-6-8-10,8': 'Cd8_2',
    '3-6-8-10,9': 'Cd8_3',
    '3-6-8-10,10': 'DC',
    '3-6-8-10,11': 'Mg_6',
    '3-6-8-10,12': 'Mg_7',
    '3-6-8-10,13': 'Pv_mac',
    '3-6-8-10,14': 'B_cell',
    '3-6-8-10,15': 'Neutrophil',
    '3-6-8-10,16': 'Inf_mac_2',
}


# name categories
immune_cat = {
    '3-6-8-10,0': 'Inf_mac',
    '3-6-8-10,1': 'Cd4',
    '3-6-8-10,2': 'Microglia',
    '3-6-8-10,3': 'Cd8',
    '3-6-8-10,4': 'Microglia',
    '3-6-8-10,5': 'Microglia',
    '3-6-8-10,6': 'Microglia',
    '3-6-8-10,7': 'Microglia',
    '3-6-8-10,8': 'Cd8',
    '3-6-8-10,9': 'Cd8',
    '3-6-8-10,10': 'DC',
    '3-6-8-10,11': 'Microglia',
    '3-6-8-10,12': 'Microglia',
    '3-6-8-10,13': 'Pv_mac',
    '3-6-8-10,14': 'B_cell',
    '3-6-8-10,15': 'Neutrophil',
    '3-6-8-10,16': 'Inf_mac',
}



#map cell types
immune.obs['cat'] = immune.obs.leiden_R.map(immune_cat)

immune.obs['cat']

##rename leiden_R names
immune.obs['leiden_R'].replace(immune_dict, inplace=True)

# save as new column in obs
immune.obs['immune_recluster'] = immune.obs['leiden_R'].copy()



# visualize
sc.pl.umap(immune, color = ['leiden_R'])

sc.tl.dendrogram(immune, groupby = "cat")

#immune_genelist = ['Rbfox3', 'Csf1r', 'Mpeg1', 'Tmem119', 'Adora3',  'Nos2', 'C3', 'Cybb', 'Cx3cr1', 'Mertk', 'Sall1', 'Cd3e', 'Cd4', 'Il2rb', 'Cxcr6', 'Cxcr3', 'Cd8a', 'Gzmb', 'Mki67', 'Ets1', 'Itgal', 'Nlrp1c-ps', 'Naip6', 'Csf3r', 'Tgfbr1', 'Tlr4', 'Mrc1', 'Siglec1', 'Cxcl9', 'Cd163', 'Axl', 'Cd19', 'Syk', 'Cd79b', 'Trem2', 'Ccr7', 'Relb', 'Il12b', 'Cd80', 'Csf2rb', 'Ifng', 'Il21', 'Tnf', 'Nr4a1', 'Csf3r', 'Tnfaip2', 'Cxcr2', 'Il1b', 'Clec7a', 'Mmp9', 'S100a9', 'Sell', 'Siglech', 'Cd47']


immune_genelist = ['Cd3e', 'Cd4', 'Cd8a', 'Gzmb', 'Mpeg1',  'Csf1r', 'Sall1', 'Tmem119', 'Il1a', 'Syk', 'Cd19', 'Cd79b', 'Xcr1', 'Ccr7', 'Relb', 'Il12b', 'Cd80', 'Mrc1', 'Cd163', 'Siglec1', 'C3', 'Il1b', 'Tnfaip2', 'Nos2', 'Chil3', 'Cybb', 'Sell', 'Mmp9', 'S100a9', 'Csf3r']
sc.pl.dotplot(immune, immune_genelist, groupby = ['cat'], dendrogram = True)


sc.pl.dotplot(immune, immune_genelist, groupby = ['immune_recluster'], dendrogram = True)




sq.pl.spatial_scatter(immune, shape=None, color=["leiden_R"], groups = 'Cd4',  wspace=0.4,)
sq.pl.spatial_scatter(immune, shape=None, color=["leiden_R"], groups = ['Cd8_1', 'Cd8_2', 'Cd8_3'],  wspace=0.4,)
sq.pl.spatial_scatter(immune, shape=None, color=["leiden_R"], groups = 'Cd8_2',  wspace=0.4,)
sq.pl.spatial_scatter(immune, shape=None, color=["leiden_R"], groups = 'Cd8_3',  wspace=0.4,)





#####################################

# T cell recluster


t_cell = adata[adata.obs['leiden'].isin(['8', '10'])]

sc.tl.pca(t_cell, svd_solver='arpack')
sc.pp.neighbors(t_cell, n_neighbors=10)
sc.tl.umap(t_cell)
sc.pl.umap(t_cell, color = ['leiden'])

sc.tl.leiden(t_cell, resolution=.5, restrict_to=('leiden', ['8', '10']))


sc.pl.umap(t_cell, color = ['leiden_R', 'Cd3e'])
          

# rank genes
sc.tl.rank_genes_groups(t_cell, 'leiden_R')
sc.pl.rank_genes_groups(t_cell, n_genes=20, sharey=False)



# name clusters

t_cell_r = {
    '8-10,0': 'Cd4_1',
    '8-10,1': 'Cd8_1',
    '8-10,2': 'Cd8_2',
    '8-10,3': 'Proliferating',
    '8-10,4': 'Cd8_3',
    '8-10,5': 'Ifng+',
    '8-10,6': 'Treg'
}


#map cell types
t_cell.obs['t_cell_r'] = t_cell.obs.leiden_R.map(t_cell_r)

t_cell.obs['t_cell_r']

##rename leiden_R names
t_cell.obs['leiden_R'].replace(t_cell_r, inplace=True)

# save as new column in obs
t_cell.obs['t_cell_r'] = t_cell.obs['t_cell_r'].copy()



# visualize
sc.pl.umap(t_cell, color = ['leiden_R'])

# t cell spatial maps
sq.pl.spatial_scatter(t_cell, shape=None, color=["t_cell_r"], groups = 'Cd4_1',  figsize = 2, wspace=0.4,)
sq.pl.spatial_scatter(t_cell, shape=None, color=["t_cell_r"], groups = 'Cd8_1',  wspace=0.4,)
sq.pl.spatial_scatter(t_cell, shape=None, color=["t_cell_r"], groups = 'Cd8_2',  wspace=0.4,)
sq.pl.spatial_scatter(t_cell, shape=None, color=["t_cell_r"], groups = 'Proliferating',  wspace=0.4,)
sq.pl.spatial_scatter(t_cell, shape=None, color=["t_cell_r"], groups = 'Ifng+',  wspace=0.4,)
sq.pl.spatial_scatter(t_cell, shape=None, color=["t_cell_r"], groups = 'Treg',  wspace=0.4,)







#####################################

#Gene scoring just immune cells


sc.pl.matrixplot(immune, 
                 var_names = ['cytokine_score', 'chemokine_score', 'gf_signaling', 'cell_activation_score', 'dam_score', 'cell_death_score'], 
                 groupby='cat', 
                 dendrogram=True,
                 standard_scale = 'var')


# Gene scoring microglia

microglia = immune[immune.obs['cat'].isin(['Microglia'])]

sc.tl.dendrogram(microglia, groupby = 'leiden_R')                  

microglia.obs['microglia_cluster'] = microglia.obs['leiden_R'].copy()

sc.pl.matrixplot(microglia, 
                 var_names = ['cytokine_score', 'chemokine_score', 'gf_signaling', 'cell_activation_score', 'dam_score', 'cell_death_score'], 
                 groupby='microglia_cluster', 
                 dendrogram=True,
                 standard_scale = 'var')


sc.tl.rank_genes_groups(microglia, 'leiden_R')
sc.pl.rank_genes_groups(microglia, n_genes=10, sharey=False)



#mask out neuronal contamination
mg_mask = microglia[~microglia.obs['leiden_R'].isin(['Mg_3'])]
sc.tl.dendrogram(mg_mask, groupby = 'microglia_cluster')                  


# recluster
sc.tl.pca(mg_mask, svd_solver = 'arpack')
sc.pp.neighbors(mg_mask, n_neighbors=10)
sc.tl.umap(mg_mask)
sc.pl.umap(mg_mask, color = ['leiden_R'])

#visualize microglia

# USE THESE FOR MG SPATIAL MAPS
sq.pl.spatial_scatter(mg_mask, shape=None, color=["immune_recluster"], groups = 'Mg_1',  wspace=0.4,)
sq.pl.spatial_scatter(mg_mask, shape=None, color=["immune_recluster"], groups = 'Mg_2',  wspace=0.4,)
sq.pl.spatial_scatter(mg_mask, shape=None, color=["immune_recluster"], groups = 'Mg_3',  wspace=0.4,)
sq.pl.spatial_scatter(mg_mask, shape=None, color=["immune_recluster"], groups = 'Mg_4',  wspace=0.4,)
sq.pl.spatial_scatter(mg_mask, shape=None, color=["immune_recluster"], groups = 'Mg_5',  wspace=0.4,)
sq.pl.spatial_scatter(mg_mask, shape=None, color=["immune_recluster"], groups = 'Mg_6',  wspace=0.4,)



sc.tl.rank_genes_groups(mg_mask, 'leiden_R')
sc.pl.rank_genes_groups(mg_mask, n_genes=20, sharey=False)


microglia_genes = ['Csf1r', 'Mpeg1', 'Sall1', 'Tmem119', 'Tgfbr2', 'Cx3cr1', 'Trem2', 'Aif1', 'Tlr12',  'Axl', 'Msr1', 'Itgax', 'Adora3', 'Naip6', 'Nlrp1c-ps', 'Tlr4', 'Tnf', 'Il1a', 'Icam1', 'Cxcl9', 'Cxcl10', 'Ccrl2', 'Cxcl16', 'Ccl2', 'Nfkbia']
sc.pl.dotplot(mg_mask, microglia_genes, groupby = 'microglia_cluster', dendrogram = True)




# gene scoring masked mg
sc.pl.matrixplot(mg_mask, 
                 var_names = ['cytokine_score', 'chemokine_score', 'gf_signaling', 'cell_activation_score', 'dam_score', 'cell_death_score'], 
                 groupby='microglia_cluster', 
                 dendrogram=True,
                 standard_scale = 'var')





#########################

# IMMUNE PLOTS

sc.pl.umap(immune, color = 'leiden_R')

sc.pl.umap(immune, color = ['Cd3e', 'Cd4', 'Cd8a'])
sc.pl.umap(immune, color = ['Mmp9'])
sc.pl.umap(immune, color = ['Sall1'])
sc.pl.umap(immune, color = ['Cd19'])
sc.pl.umap(immune, color = ['Csf1r'])
sc.pl.umap(immune, color = ['Tnfaip2'])
sc.pl.umap(immune, color = ['Nos2'])
sc.pl.umap(immune, color = ['Mrc1'])
sc.pl.umap(immune, color = ['Mki67'])
sc.pl.umap(immune, color = ['Syk'])
sc.pl.umap(immune, color = ['Mpeg1'])


sc.pl.umap(immune, color = ['Csf1r', 'Cd19', 'Xcr1', 'Ccr7'])


immune_genelist = ['Cd3e', 'Cd4', 'Cd8a', 'Gzmb', 'Mpeg1',  'Csf1r', 'Sall1', 'Tmem119', 'Il1a', 'Syk', 'Cd19', 'Cd79b', 'Xcr1', 'Ccr7', 'Relb', 'Il12b', 'Cd80', 'Mrc1', 'Cd163', 'Siglec1', 'C3', 'Il1b', 'Tnfaip2', 'Nos2', 'Chil3', 'Cybb', 'Sell', 'Mmp9', 'S100a9', 'Csf3r']
sc.pl.tracksplot(immune, var_names = immune_genelist, groupby = 'cat',  show_gene_labels=True, dendrogram = True)



cell_death2 = ['Fasl', 'Casp1', 'Casp3', 'Casp4',  'Casp7', 'Casp8', 'Nlrp3', 'Nlrc3', 'Nlrp1a', 'Nlrp1b', 'Nlrp12', 'Naip2', 'Nlrp1c-ps', 'Gsdmd', 'Pycard', 'Tnfrsf1a', 'Tnfrsf1b', 'Ripk1', 'Ripk3', 'Fas',  'Tnfrsf21', 'Apaf1', 'Cflar', 'Birc2', 'Birc3', 'Traf2',  'Nfkbia', 'Ticam1', 'Zbp1', 'Bak1', 'Bcl2l11', 'Pmaip1', 'Blk',  'S100a9', 'Tnf']
sc.pl.dotplot(immune, var_names = cell_death2, groupby = 'cat', dendrogram = True)



sc.pl.dotplot(immune, immune_genelist, groupby = ['immune_recluster'], dendrogram = True)






##########################################

# reclustering immune cells without contamination in mg_3

immune2 = immune[~immune.obs['leiden_R'].isin(['Mg_3'])]

sc.tl.pca(immune2, svd_solver='arpack')
sc.pp.neighbors(immune2, n_neighbors=10)
sc.tl.umap(immune2)
sc.pl.umap(immune2, color = ['leiden'])

sc.tl.leiden(immune2, resolution=.8, restrict_to=('leiden', ['3', '6', '8', '10']))


          
sc.tl.rank_genes_groups(immune2, 'leiden_R')
sc.pl.rank_genes_groups(immune2, n_genes=20, sharey=False)


# Rename clusters

immune_dict2 = {
    '3-6-8-10,0': 'Inf_mac_1',
    '3-6-8-10,1': 'Mg_1',
    '3-6-8-10,2': 'Cd4',
    '3-6-8-10,3': 'Cd8_1',
    '3-6-8-10,4': 'Mg_2',
    '3-6-8-10,5': 'Mg_3',
    '3-6-8-10,6': 'Mg_4',
    '3-6-8-10,7': 'Mg_5',
    '3-6-8-10,8': 'Prolif_T_cell',
    '3-6-8-10,9': 'Mg_6',
    '3-6-8-10,10': 'Cd8_2',
    '3-6-8-10,11': 'Pv_mac',
    '3-6-8-10,12': 'DC',
    '3-6-8-10,13': 'B_cell',
    '3-6-8-10,14': 'Neutrophil',
    '3-6-8-10,15': 'Inf_mac',
}


# name categories
immune_cat2 = {
    '3-6-8-10,0': 'Inf_mac',
    '3-6-8-10,1': 'Microglia',
    '3-6-8-10,2': 'Cd4',
    '3-6-8-10,3': 'Cd8',
    '3-6-8-10,4': 'Microglia',
    '3-6-8-10,5': 'Microglia',
    '3-6-8-10,6': 'Microglia',
    '3-6-8-10,7': 'Microglia',
    '3-6-8-10,8': 'Prolif_T_cell',
    '3-6-8-10,9': 'Microglia',
    '3-6-8-10,10': 'Cd8',
    '3-6-8-10,11': 'Pv_mac',
    '3-6-8-10,12': 'DC',
    '3-6-8-10,13': 'B_cell',
    '3-6-8-10,14': 'Neutrophil',
    '3-6-8-10,15': 'Inf_mac'
    }



#map cell types
immune2.obs['cat'] = immune2.obs.leiden_R.map(immune_cat2)

immune2.obs['cat'] # check before running the next line

##rename leiden_R names
immune2.obs['leiden_R'].replace(immune_dict2, inplace=True)

# save as new column in obs
immune2.obs['immune_recluster'] = immune2.obs['leiden_R'].copy()



# visualize
sc.pl.umap(immune2, color = ['cat'])

sc.tl.dendrogram(immune, groupby = "immune_recluster")

immune_genelist = ['Cd3e', 'Cd4', 'Cd8a', 'Gzmb', 'Mpeg1',  'Csf1r', 'Sall1', 'Tmem119', 'Il1a', 'Syk', 'Cd19', 'Cd79b', 'Xcr1', 'Ccr7', 'Relb', 'Il12b', 'Cd80', 'Mrc1', 'Cd163', 'Siglec1', 'C3', 'Il1b', 'Tnfaip2', 'Nos2', 'Chil3', 'Cybb', 'Sell', 'Mmp9', 'S100a9', 'Csf3r']
sc.pl.dotplot(immune, immune_genelist, groupby = ['immune_recluster'], dendrogram = True)





#########


sc.pl.heatmap(immune, cell_death, groupby = ['cat'], dendrogram = True)

##########################################

# IMMUNE SPATIAL MAPS

## group names
Microglia = ['Mg_1', 'Mg_2', 'Mg_3', 'Mg_4', 'Mg_5', 'Mg_6', 'Mg_7']
Inf_mac = ['Inf_mac_1', 'Inf_mac_2']
Cd8 = ['Cd8_1', 'Cd8_2']

# Spatial maps
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = ['B_cell', 'Pv_mac', 'DC'],  wspace=0.4,)

sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'B_cell',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'Neutrophil',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'Pv_mac',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'DC',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'Inf_mac',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'Cd8',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'Cd4',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["cat"], groups = 'Microglia',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["immune_recluster"], groups = 'DC',  wspace=0.4,)


## microglia
### when reclustered with immune cells heterogenety change
sq.pl.spatial_scatter(immune2, shape=None, color=["immune_recluster"], groups = 'Mg_1',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["immune_recluster"], groups = 'Mg_2',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["immune_recluster"], groups = 'Mg_3',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["immune_recluster"], groups = 'Mg_4',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["immune_recluster"], groups = 'Mg_5',  wspace=0.4,)
sq.pl.spatial_scatter(immune2, shape=None, color=["immune_recluster"], groups = 'Mg_6',  wspace=0.4,)

#########


# NEIGHBORHOOD ENRICHMENT

sq.gr.spatial_neighbors(immune2)
sq.gr.nhood_enrichment(immune2, cluster_key="cat")
sq.pl.nhood_enrichment(immune2, cluster_key="cat", cmap = 'inferno')


sq.gr.co_occurrence(immune2, cluster_key="cat")
sq.pl.co_occurrence(immune2, cluster_key="cat", clusters = 'Pv_mac')

sq.gr.interaction_matrix(immune2, cluster_key = 'cat')
sq.pl.interaction_matrix(immune2, cluster_key = 'cat')
#########



# microglia scoring fig

#microglia2 = immune2[immune2.obs['cat'].isin(['Microglia'])]

#sc.tl.dendrogram(microglia2, groupby = 'immune_recluster')                  

#microglia2.obs['microglia_cluster'] = microglia2.obs['leiden_R'].copy()

#sc.pl.matrixplot(microglia2, 
#                 var_names = ['cytokine_score', 'chemokine_score', 'gf_signaling', 'cell_activation_score', 'dam_score', 'cell_death_score'], 
#                 groupby='microglia_cluster', 
 #                dendrogram=True,
 #                standard_scale = 'var')

##########################################

# DRAFT Macrophage recluster 

macrophage = adata[adata.obs['leiden'].isin(['3', '6'])]

sc.tl.pca(macrophage, svd_solver='arpack')
sc.pp.neighbors(macrophage, n_neighbors=15)
sc.tl.umap(macrophage)
sc.pl.umap(macrophage, color = ['leiden'])

sc.tl.leiden(macrophage, resolution=.8, restrict_to=('leiden', ['3', '6']))


sc.pl.umap(macrophage, color = ['leiden_R', 'Csf1r', 'Mmp9', 'S100a9', 'Mrc1', 'Xcr1', 'Itgax', 'Cd19'])



# rank genes
sc.tl.rank_genes_groups(macrophage, 'leiden_R')
sc.pl.rank_genes_groups(macrophage, n_genes=10, sharey=False)


# Rename clusters

macrophage_dict = {
    '3-6,0': 'Mg_1',
    '3-6,1': 'Inf_mac_1',
    '3-6,2': 'Mg_2',
    '3-6,3': 'Mg_3',
    '3-6,4': 'Pv_mac',
    '3-6,5': 'B_cell',
    '3-6,6': 'Unknown',
    '3-6,7': 'Neutrophil',
    '3-6,8': 'Mixed'
}

macrophage.obs['leiden_R'].replace(macrophage_dict, inplace=True)


# save as new column in obs
macrophage.obs['mac_recluster'] = macrophage.obs['leiden_R'].copy()




#subset out contaminants

macrophage_filtered = macrophage[macrophage.obs['leiden_R'].isin(['Mg_1', 'Mg_2', 'Mg_3', 'Pv_mac', 'Inf_mac_1'])]


sc.pl.umap(macrophage_filtered, color = ['mac_recluster', 'P2ry12', 'Cd163'])



macrophage = macrophage[~macrophage.obs['leiden_R'].isin(['3-6,3', '3-6,6', '3-6,8', '3-6,9', '3-6,10'])]



sc.pl.umap(macrophage, color = ['leiden', 'leiden_R', 'Il1a', 'P2ry12', 'Tnfaip2', 'Csf2rb2', 'Nlrp1c-ps','C3', 'Tnfaip2', 'Nos2', 'Clec7a',  'Apoc2', 'Csf1r', 'Axl'])
sc.tl.rank_genes_groups(macrophage, 'leiden_R')
sc.pl.rank_genes_groups(macrophage, n_genes=14, sharey=False)
sc.tl.dendrogram(macrophage, groupby = 'leiden_R')





mac_types= {
"3-6,0": "Infiltrating macrophage 1",
"3-6,1": "Microglia 1",
"3-6,2": "Microglia 2",
"3-6,4": "Microglia 3",
"3-6,5": "Infiltrating macrophage 2",
"3-6,7": "Perivascular macrophage"
}

macrophage.obs['mac_types'] = macrophage.obs.leiden.map(mac_types)


#mac_markers = ['Axl', 'Mertk', 'Xcr1', 'Cd86', 'Cd79b', 'Il1b', 'Il1a', 'Stat3', 'Csf2rb2','Il1b', 'S100a9', 'Sell', 'Stat1', 'Csf3r', 'Nlrp3', 'Syk', 'Tlr11', 'Itgax', 'Chil3', 'Tmem119', 'Hpgds', 'P2ry12', 'Il21r', 'Csf1r', 'Csf2ra', 'Adora3', 'Relb', 'Cx3cr1', 'Sall1', 'Tnfaip2', 'Mpeg1', 'Siglech', 'C3', 'Cybb', 'Tnfaip2', 'Cxcl9', 'Cxcl10', 'Clec7a', 'Cxcl16', 'Tnf', 'Nos2', 'Ccrl2', 'Adora3', 'Naip6', 'Tlr4', 'Trem2', "Nlrp12", 'Nlrp1c-ps', 'Nlrp1b', 'Itga4', 'Ccr7']
#sc.pl.dotplot(macrophage, mac_markers, groupby = 'leiden_R', dendrogram = True)


mac2 = ['Mpeg1', 'Csf1r', 'Tmem119', 'Sall1', 'Trem2', 'Naip6', 'Adora3', 'Nlrp1c-ps',  'Tnfaip2', 'C3', 'Cybb', 'Icam1', 'Cxcl10','Relb','Il1b', 'Ccl2', 'Il1a', 'Tnf',  'Clec7a',  'Itga4',   'Nos2', 'Ccr2',   'Mrc1', 'Siglec1', 'Cd163' ]
sc.pl.dotplot(macrophage, mac2, groupby = 'leiden_R', dendrogram = True)

sc.pl.umap(macrophage, color = ['leiden_R', 'Il1a'])
sc.pl.umap(macrophage, color = ['leiden_R', 'Mrc1', 'Tmem119', 'Nos2'])



#############################

mg = ['3-6,1', '3-6,2', '3-6,4']
sq.pl.spatial_scatter(macrophage, shape=None, color=["leiden_R"], groups = mg,  wspace=0.4,)



mac = ['3-6,5', '3-6,0']
sq.pl.spatial_scatter(macrophage, shape=None,color=["leiden_R"], groups = mac, wspace=0.4,)


pv = ['3-6,7']
sq.pl.spatial_scatter(macrophage, shape=None,color=['leiden_R'], groups = pv, wspace=0.4,)
 


sq.pl.spatial_scatter(macrophage, shape = None, color = ['Tmem119', 'Tnfaip2', 'Cxcl9', 'C3'], outline = True)


sq.gr.spatial_neighbors(macrophage, key_added = 'spatial')
sq.gr.nhood_enrichment(macrophage, cluster_key = 'leiden')
sq.pl.nhood_enrichment(macrophage, cluster_key = 'leiden')

sq.gr.interaction_matrix(macrophage, cluster_key = 'leiden_R')
sq.pl.interaction_matrix(macrophage, cluster_key = 'leiden_R')
 
sq.gr.ripley(macrophage, mode = 'L', cluster_key = 'leiden')   
sq.pl.ripley(macrophage,  cluster_key = 'leiden')

sq.gr.co_occurrence(macrophage, cluster_key = 'leiden_R')
sq.pl.co_occurrence(macrophage, 'leiden_R')

    


cell_counts = macrophage.obs['leiden_R'].value_counts()
cell_counts = pd.DataFrame(cell_counts)
cell_counts['leiden_R'] = cell_counts.index
sns.barplot(cell_counts, x = 'count', y = 'Cell_type')



####

sq.pl.spatial_scatter(adata, color=["Ifng", "cluster"])

