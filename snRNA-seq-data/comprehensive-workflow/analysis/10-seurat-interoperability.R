library(Seurat)

################################################################################################################################

# FULL DATASET, ALL GENES

################################################################################################################################

# SET UP DIRECTORIES

base_dir <- '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/csv/full'

## Define paths for csv exports of original AnnData
counts_file <- file.path(base_dir, 'data', 'counts-matrix.csv')
log1p_file <- file.path(base_dir, 'data', 'log1p-matrix.csv')
normalized_file <- file.path(base_dir, 'data', 'normalized-matrix.csv')
cell_metadata_file <- file.path(base_dir, 'metadata', 'cell-metadata.csv')
feature_metadata_file <- file.path(base_dir, 'metadata', 'feature-metadata.csv')
pca_file <- file.path(base_dir, 'reductions', 'pca.csv')
latent_rep_file <- file.path(base_dir, 'reductions', 'latent_representation.csv')
umap_file <- file.path(base_dir, 'reductions', 'umap.csv')

################################################################################################################################

# MANUALLY ASSEMBLE SEURAT OBJECT (FULL)

## Load counts matrix
counts_matrix <- read.csv(counts_file, row.names = 1)
counts_matrix <- t(as.matrix(counts_matrix)) # transpose

## Load normalized and log1p count matrices
normalized_matrix <- read.csv(normalized_file, row.names = 1)
normalized_matrix <- t(as.matrix(normalized_matrix))

log1p_matrix <- read.csv(log1p_file, row.names = 1)
log1p_matrix <- t(as.matrix(log1p_matrix))

## Create SeuratObject
seurat_object <- CreateSeuratObject(counts = counts_matrix, project = "2024_snRNA-seq_TBI")

seurat_object[["log1p"]] <- CreateAssayObject(counts = log1p_matrix)
seurat_object[["normalized"]] <- CreateAssayObject(counts = normalized_matrix)

## Load and add cell metadata
cell_metadata <- read.csv(cell_metadata_file, row.names = 1)
seurat_object <- AddMetaData(seurat_object, metadata = cell_metadata)

## Load and add feature (gene) metadata
feature_metadata <- read.csv(feature_metadata_file, row.names = 1)

common_features <- intersect(rownames(seurat_object@assays$RNA), rownames(feature_metadata))
if (length(common_features) == 0) {
  stop("No feature overlap between counts matrix and feature metadata.")
}
feature_metadata <- feature_metadata[common_features, ]

for (colname in colnames(feature_metadata)) {
  seurat_object <- AddMetaData(seurat_object, metadata = feature_metadata[[colname]], col.name = colname)
}

## Load and add PCA coordinates
pca_coords <- read.csv(pca_file, check.names = FALSE, row.names = 1)
colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))
pca_coords <- as.matrix(pca_coords)
seurat_object[['pca']] <- CreateDimReducObject(embeddings = pca_coords, key = "PC_", assay = DefaultAssay(seurat_object))

## Load and add scVI latent representation
latent_rep <- read.csv(latent_rep_file, check.names = FALSE, row.names = 1)
colnames(latent_rep) <- paste0("latent_", 1:ncol(latent_rep))
latent_rep <- as.matrix(latent_rep)
seurat_object[['latent']] <- CreateDimReducObject(embeddings = latent_rep, key = "latent_", assay = DefaultAssay(seurat_object))

## Load and add UMAP coordinates
umap_coords <- read.csv(umap_file, row.names = 1)
colnames(umap_coords) <- paste0("umap_", 1:ncol(umap_coords))
umap_coords <- as.matrix(umap_coords)
seurat_object[['umap']] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = DefaultAssay(seurat_object))

################################################################################################################################

# DATA VIZ (FULL)

## UMAP of annotations
DimPlot(seurat_object, reduction = 'umap', label = TRUE, pt.size = 0.5, group.by = 'cell_type')
DimPlot(seurat_object, reduction = 'umap', label = TRUE, pt.size = 0.5, group.by = 'cluster')
DimPlot(seurat_object, reduction = 'umap', label = TRUE, pt.size = 0.5, group.by = 'leiden')

################################################################################################################################

# EXPORT DATASET (FULL)

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/seurat'

saveRDS(seurat_object, file = file.path(save_dir, "2024_tbi-snrna-seq-cleaned-full.rds"))


################################################################################################################################
################################################################################################################################

# HIGHLY-VARIABLE GENES SUBSET

################################################################################################################################

# SET UP DIRECTORIES

base_dir <- '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/csv/hvg'

## Define paths for csv exports of original AnnData
counts_file <- file.path(base_dir, 'data', 'counts-matrix.csv')
log1p_file <- file.path(base_dir, 'data', 'log1p-matrix.csv')
normalized_file <- file.path(base_dir, 'data', 'normalized-matrix.csv')
cell_metadata_file <- file.path(base_dir, 'metadata', 'cell-metadata.csv')
feature_metadata_file <- file.path(base_dir, 'metadata', 'feature-metadata.csv')
pca_file <- file.path(base_dir, 'reductions', 'pca.csv')
latent_rep_file <- file.path(base_dir, 'reductions', 'latent_representation.csv')
umap_file <- file.path(base_dir, 'reductions', 'umap.csv')

################################################################################################################################

# MANUALLY ASSEMBLE SEURAT OBJECT (HVG)

## Load counts matrix
counts_matrix <- read.csv(counts_file, row.names = 1)
counts_matrix <- t(as.matrix(counts_matrix))

## Load normalized and log1p count matrices
normalized_matrix <- read.csv(normalized_file, row.names = 1)
normalized_matrix <- t(as.matrix(normalized_matrix))

log1p_matrix <- read.csv(log1p_file, row.names = 1)
log1p_matrix <- t(as.matrix(log1p_matrix))


## Create SeuratObject
seurat_object <- CreateSeuratObject(counts = counts_matrix, project = "2024_snRNA-seq_TBI")

seurat_object[["log1p"]] <- CreateAssayObject(counts = log1p_matrix)
seurat_object[["normalized"]] <- CreateAssayObject(counts = normalized_matrix)

## Load and add cell metadata
cell_metadata <- read.csv(cell_metadata_file, row.names = 1)
seurat_object <- AddMetaData(seurat_object, metadata = cell_metadata)

## Load and add feature (gene) metadata
feature_metadata <- read.csv(feature_metadata_file, row.names = 1)

common_features <- intersect(rownames(seurat_object@assays$RNA), rownames(feature_metadata))
if (length(common_features) == 0) {
  stop("No feature overlap between counts matrix and feature metadata.")
}
feature_metadata <- feature_metadata[common_features, ]

for (colname in colnames(feature_metadata)) {
  seurat_object <- AddMetaData(seurat_object, metadata = feature_metadata[[colname]], col.name = colname)
}


## Load and add PCA coordinates
pca_coords <- read.csv(pca_file, check.names = FALSE, row.names = 1)
colnames(pca_coords) <- paste0("PC_", 1:ncol(pca_coords))
pca_coords <- as.matrix(pca_coords)
seurat_object[['pca']] <- CreateDimReducObject(embeddings = pca_coords, key = "PC_", assay = DefaultAssay(seurat_object))

## Load and add scVI latent representation
latent_rep <- read.csv(latent_rep_file, check.names = FALSE, row.names = 1)
colnames(latent_rep) <- paste0("latent_", 1:ncol(latent_rep))
latent_rep <- as.matrix(latent_rep)
seurat_object[['latent']] <- CreateDimReducObject(embeddings = latent_rep, key = "latent_", assay = DefaultAssay(seurat_object))

## Load and add UMAP coordinates
umap_coords <- read.csv(umap_file, row.names = 1)
colnames(umap_coords) <- paste0("umap_", 1:ncol(umap_coords))
umap_coords <- as.matrix(umap_coords)
seurat_object[['umap']] <- CreateDimReducObject(embeddings = umap_coords, key = "UMAP_", assay = DefaultAssay(seurat_object))

################################################################################################################################

# DATA VIZ (HVG)

## Plot UMAP
DimPlot(seurat_object, reduction = 'umap', label = TRUE, pt.size = 0.5, group.by = 'cell_type')
DimPlot(seurat_object, reduction = 'umap', label = TRUE, pt.size = 0.5, group.by = 'cluster')
DimPlot(seurat_object, reduction = 'umap', label = TRUE, pt.size = 0.5, group.by = 'leiden')

################################################################################################################################

# EXPORT (HVG)

save_dir = '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq/analysis/1/export/seurat'

saveRDS(seurat_object, file = file.path(save_dir, "2024_tbi-snrna-seq-cleaned-hvg.rds"))
