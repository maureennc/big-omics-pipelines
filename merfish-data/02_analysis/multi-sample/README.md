# MERFISH Analysis

## Overview
This repository contains scripts for processing and analyzing MERFISH (Multiplexed Error-Robust Fluorescence In Situ Hybridization) data, focusing on spatially-resolved gene expression patterns and immune responses in various brain regions.

## Scripts

### `01_pp-qc-concat.py`
- Import and process datasets from multiple MERFISH experiments (3 samples).
- Apply quality control thresholds, score doublets using Scrublet, and filter out low-quality cells.
- Concatenate, normalize, and save the processed datasets in H5AD format for further analysis.

### `02_pp-polar-coordinates.py`
- Process concatenated MERFISH data by filtering out doublets using predefined thresholds.
- Calculate spatial centroids for each cell and convert these coordinates into polar coordinates (radius and angle).
- Evaluate multicollinearity among several covariates through correlation analysis, visualized with a heatmap.
- Export the modified AnnData object containing updated spatial information and filtering results for further analysis.

### `03_scvi-training-annotation.py`
- Process and analyze MERFISH data using the scVI model to infer latent representations and cluster cell populations.
- Import processed AnnData objects and set random seeds for reproducibility.
- Set up the scVI model with appropriate covariates, including categorical variables (sample, condition) and continuous variables such as total counts, doublet scores, and polar coordinates (centroid theta and radius).
- Train the model and evaluate its performance by extracting and visualizing the training history, including metrics such as ELBO and training loss over epochs.
- Perform initial clustering using the latent representation and visualize results in UMAP space.
- Identify marker genes for each cluster and validate clustering quality by analyzing doublet scores.
- Annotate cell types and clusters based on clustering results, and export the annotated AnnData object in H5AD format for further analysis.

### `04_vizualizer-split.py`
- After data cleaning, clustering, and cell type annotation, split the concatenated AnnData object into individual datasets for each sample, as the `Vizualizer` can only import one experiment sample at a time. Remove barcode suffixes to ensure compatibility, matching the original experiment barcodes.
- Use the `Vizualizer` software to draw regions of interest (ROIs) around inflammatory foci. Export polygon geometries for these immune hotspots (x-y coordinates for vertices) for further downstream analysis in Python.
- Identify regions of interest based on polygon coordinates defined by exported geometries (.csv). Calculate the geometric centers of cell coordinates to determine if each cell falls within the specified ROIs.
- Annotate cells within specific inflammatory cluster geometries and naive brain regions. Finally, re-concatenate the annotated datasets and export the updated AnnData objects in H5AD format for further analysis.

### `05_reintegration.py`
- Train two integration models with different combinations of three MERFISH samples: E003 (infected), E007 (infected), and E008 (naive).
- The first model integrates all three samples (essentially a concatenation) to capture the overall cellular landscape.
- The second model focuses solely on the two infected samples to potentially enhance model performance.
- This step is necessary due to previous splitting of the samples for Vizualizer compatibility. Although the samples could have been concatenated directly, splitting was performed to adhere to the limitations of the Vizgen Vizualizer software, which only allows for one experiment sample at a time.

### `06_spatial-de-clusters.py`
- Perform spatially-resolved differential expression analysis, focusing on immune cell types within spatial regions of interest (ROIs) defined by the `inflammatory_foci` annotation.
- Use two AnnData objects: `adata_concat`, which encompasses all samples for broad comparisons, and `adata_infected`, which targets differential expression specifically within infected samples.
- Conduct differential expression analyses using the Wilcoxon rank-sum test, providing insights into gene activity across specified conditions.
- Focus on immune cell types and inflammatory foci, crucial for understanding the spatial dynamics of immune responses in tissues.
- Generate visualizations to depict gene expression differences, highlighting significant findings in the context of spatial biology.

### `07_spatial-de-naive.py`
- Utilize the naive sample of MERFISH data to perform spatially-resolved differential expression analysis across different brain region annotations.
- Assess immune response profiles in various regions, revealing strong priming of the choroid plexus for immune responses.

- **Brain Regions**:
  - Cortex
  - Hippocampus
  - Striatum
  - Choroid Plexus
  
- Import the naive sample annotated with the aforementioned regional brain regions for analysis.
- Define several gene lists, including cell death genes, cytokines, and other relevant gene sets for evaluating immune responses.
- Employ the Wilcoxon rank-sum test to identify differentially expressed genes across the specified brain regions.
- Visualize results using matrix plots and bar plots to showcase the significant gene counts per region.
- Highlight strong immune activation in the choroid plexus, emphasizing its role in the immune response.

This script provides insights into the spatial context of immune cell activity in the brain, which is critical for understanding neuroimmune interactions in health and disease.

---

## Data Availability
A link to the dataset in GEO will be added soon.

---
