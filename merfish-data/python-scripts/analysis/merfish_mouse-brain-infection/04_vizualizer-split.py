#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  3 14:06:17 2024

@author: maureen
"""

import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.sparse import csr_matrix
from copy import deepcopy
from shapely.geometry import Point, Polygon
from shapely.wkt import loads

# IMPORT DATA

data_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata = sc.read_h5ad(os.path.join(data_dir, 'scvi-trained_A-annotated.h5ad'))

###############################################################################

# SPLIT SAMPLES FOR INTERACTIVE VIEWING

## Create subsets
E003 = adata[adata.obs['sample'] == 'E003-infected'].copy()
E007 = adata[adata.obs['sample'] == 'E007-infected'].copy()
E008 = adata[adata.obs['sample'] == 'E008-naive'].copy()

## Remove barcode suffix for compatibility with Vizualizer
E003.obs.index = E003.obs.index.str.replace("-E003-infected", "", regex=False)
E007.obs.index = E007.obs.index.str.replace("-E007-infected", "", regex=False)
E008.obs.index = E008.obs.index.str.replace("-E008-infected", "", regex=False)

###############################################################################

# EXPORT

save_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/vizualizer/annotated_h5ad'

E003.X = csr_matrix(E003.X)
E007.X = csr_matrix(E007.X)
E008.X = csr_matrix(E008.X)

#E003.write_h5ad(os.path.join(save_dir, 'E003-trained_A-annotated.hdf5'))
#E007.write_h5ad(os.path.join(save_dir, 'E007-trained_A-annotated.hdf5'))
#E008.write_h5ad(os.path.join(save_dir, 'E008-trained_A-annotated.hdf5'))

###############################################################################

# CREATE REGIONS OF INTEREST USING VIZGEN VIZUALIZER SOFTWARE

## Use polygon tool in GUI to draw custom ROI and create polygon geometry csv files for compatibility with AnnData

###############################################################################

# IMPORT GEOMETRY CSV FILES

coord_dir = "/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/vizualizer/export/roi-coordinates"

## E003-infected
E003_roi = pd.read_csv(os.path.join(coord_dir, 'E003-polygon-coordinates.csv'))

E003_crop = E003_roi[E003_roi['group'] == 'crop'].copy()
E003_clusters = E003_roi[E003_roi['group'] == 'clusters'].copy()
                
## E007-infectd
E007_roi = pd.read_csv(os.path.join(coord_dir, 'E007-polygon-coordinates.csv'))

E007_crop = E007_roi[E007_roi['group'] == 'crop'].copy()
E007_clusters = E007_roi[E007_roi['group'] == 'clusters'].copy()

## E008-naive
E008_regions = pd.read_csv(os.path.join(coord_dir, 'E008-polygon-coordinates.csv'))

###############################################################################

# CROP FOLDED TISSUE

## E003
E003_crop['geometry'] = E003_crop['geometry'].apply(loads)

E003.obs['cropped'] = False

for _, row in E003_crop.iterrows():
    polygon = row['geometry']
    is_cropped = [polygon.contains(Point(x, y)) for x, y in E003.obsm['spatial']]  # Create a boolean mask for cells within the current polygon
    
    E003.obs['cropped'] |= is_cropped ### Now, E003.obs['cropped'] is True for cells inside any of the polygons, and False otherwise

E003 = E003[E003.obs['cropped'] == True].copy()



## E007
E007_crop['geometry'] = E007_crop['geometry'].apply(loads)

E007.obs['cropped'] = False

for _, row in E007_crop.iterrows():
    polygon = row['geometry']
    is_cropped = [polygon.contains(Point(x, y)) for x, y in E007.obsm['spatial']]  
    
    E007.obs['cropped'] |= is_cropped 

E007 = E007[E007.obs['cropped'] == True].copy()


## E008
E008.obs['cropped'] = False


###############################################################################

# ANNOTATE CELLS WITHIN INFLAMMATORY CLUSTER GEOMETRIES

## E003
E003_clusters['geometry'] = E003_clusters['geometry'].apply(loads)

E003.obs['inflammatory_foci'] = False

for _, row in E003_clusters.iterrows():
    polygon = row['geometry']
    is_inflammatory_foci = [polygon.contains(Point(x, y)) for x, y in E003.obsm['spatial']]  
    
    E003.obs['inflammatory_foci'] |= is_inflammatory_foci 


## E007
E007_clusters['geometry'] = E007_clusters['geometry'].apply(loads)

E007.obs['inflammatory_foci'] = False

for _, row in E007_clusters.iterrows():
    polygon = row['geometry']
    is_inflammatory_foci = [polygon.contains(Point(x, y)) for x, y in E007.obsm['spatial']] 
    
    E007.obs['inflammatory_foci'] |= is_inflammatory_foci 


## E008
E008.obs['inflammatory_foci'] = False

###############################################################################

# ANNOTATE CELLS WITHIN NAIVE BRAIN REGIONS

E008_regions['geometry'] = E008_regions['geometry'].apply(loads)

E008.obs['region_of_interest'] = ""

for _, row in E008_regions.iterrows():
    polygon = row['geometry']
    entity_id = row['EntityID'].split('-')[0]  # Strip the "-#" suffix from the EntityID
    is_in_region = [polygon.contains(Point(x, y)) for x, y in E008.obsm['spatial']]
    
    E008.obs.loc[is_in_region, 'region_of_interest'] = entity_id

E008.obs['region_of_interest'] = E008.obs['region_of_interest'].replace("", np.nan)


E008.obs.region_of_interest.value_counts(dropna=False)

###############################################################################

# RE-CONCATENATE ANNDATA FROM SPLIT ANNOTATED SAMPLES

## Address suffices first
E003.obs.index = E003.obs.index.str.replace('-E003-infected', '')
E007.obs.index = E007.obs.index.str.replace('-E007-infected', '')
E008.obs.index = E008.obs.index.str.replace('-E008-naive', '')

print(E003.obs.head())
print(E007.obs.head())
print(E008.obs.head())


## All samples
adata = E003.concatenate(E007, E008, batch_key='batch', batch_categories=['E003', 'E007', 'E008'] )
adata.obs.head()


## Infected samples
infected = E003.concatenate(E007, batch_key='sample', batch_categories=['E003', 'E007'] )
infected.obs.head()

###############################################################################

# EXPORT

save_dir = '/Users/maureen/Documents/experiments/merfish/analysis/E003-E007-E008/h5ad/processed-data'

adata.write_h5ad(os.path.join(save_dir, 'concat-model_A-roi-annotated.h5ad'))
infected.write_h5ad(os.path.join(save_dir, 'infected-model_A-roi-annotated.h5ad'))
E008.write_h5ad(os.path.join(save_dir, 'naive-model_A-roi-annotated.h5ad'))


## test import
#test = sc.read_h5ad(os.path.join(save_dir, 'naive-model_A-roi-annotated.h5ad'))

###############################################################################


