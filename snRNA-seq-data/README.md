# snRNA-seq Data Analysis Workflows

This repository contains snRNA-seq workflows analyzing two Alzheimer’s disease mouse model datasets. 

---

## Folder Structure

- **[basic-workflow/](basic-workflow/README.md)**: 
  - Streamlined pipeline for the **genetic knockout** model, focused on rapid data processing and efficient analysis to glean major biological insights (2 groups).

- **[advanced-workflow/](advanced-workflow/README.md)**: 
  - Comprehensive workflow for the **TBI + pharmacological rescue** model, featuring advanced QC, scVI modeling, and differential expression testing (6 groups).

---

## Basic Workflow

### Key Steps:
- **QC**: Basic filtering of low-quality cells.
- **Normalization**: Log-transformation and scaling.
- **Clustering**: Cell type clustering and DE.
- **Visualization**: UMAP projections to explore cell populations.

See [basic-workflow/README.md](basic-workflow/README.md)

---

## Advanced Workflow

### Key Steps:
- **Advanced QC**: Uses quantile-based thresholds and Scrublet for rigorous doublet detection.
- **scVI Modeling**: Trains multiple scVI models and integrates highly variable genes with full-genome data for downstream differential expression analysis.
- **Differential Expression**: Conducts Wald test-based differential expression testing between treatment groups, generating detailed heatmaps and visualizations.
- **Seurat Interoperability**: Prepares data layers for easy integration with Seurat in R, facilitating cross-environment analysis.

See the full workflow in the [advanced-workflow/README.md](advanced-workflow/README.md).

---

For more details, refer to the **README.md** files within each folder.
