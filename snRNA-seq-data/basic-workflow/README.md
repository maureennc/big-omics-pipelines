---

# Basic snRNA-seq Analysis Workflow

## Background

### Experimental Design
This basic workflow uses snRNA-seq data from a mouse brain experiment involving two groups. The data comes from three separate experiments focused on different types of genetic knockouts (WT, global KO, and two conditional cre KO) for a gene of interest in an Alzheimer's mouse model.

---

## Scripts

This section focuses on Python scripts since it assumes Cell Ranger has already been run on your samples (if you are using Azenta for sequencing). If you are starting from `.fastq` files instaed of `cellranger` outs, see the following scripts:

- [scRNA-seq-cellranger.sh](../../scRNA-seq-data/standard-workflow/01_cellranger.sh)
- [snRNA-seq-cellranger.sh](../../snRNA-seq-data/advanced-workflow/01_cellranger.sh)

---

### `01_pre-processing.py`
Processes WT and global KO snRNA-seq data, applying quality control thresholds, Scrublet doublet detection, and concatenation of datasets. Selects highly variable genes and genes of interest, performs clustering and UMAP visualization, and exports the processed data.


### `02_integration.py`
Integrates experimental and control snRNA-seq data using scVI, performing dimensionality reduction, clustering, and marker gene analysis. Trains two SCVI models with iterative filtering, refines cluster annotations, evaluates model performance, and exports the processed data.


### `03_differential-expression.py`
Conducts Wald test-based differential expression analysis with median-based normalization and log1p transformation at both cell type and cluster levels. Filters results based on significance criteria and identifies upregulated genes.


### `04_de-summary.py`
Visualizes the number of differentially expressed genes (DEGs) across cell types and clusters. Generates bar plots to show DEG distribution based on filtered results (q-value < 0.05, |log2FC| ≥ 0.5).


### `05_volcano-plots.py`
Generates volcano plots for multiple cell types using DE results. Highlights specific genes of interest and saves plots in both PDF and PNG formats.


### `06_export-anndata-components.py`
Exports AnnData components into CSV files for interoperability with Seurat in R.


### `07_create-rds-seurat.R`
Constructs a Seurat object in R using the exported CSV files from the AnnData object.

---

## Data Availability
A link to the dataset in GEO will be added soon.

---
