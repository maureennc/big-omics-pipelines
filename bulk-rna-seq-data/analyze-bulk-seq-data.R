# ANALYZE BULK RNA-SEQ DATA

#######################################################################################################

# PACKAGES

library(biomaRt)
library(tximport)
library(DESeq2)
library(tibble)
library(topGO)
require(topGO)
require(org.Mm.eg.db)

library(ggplot2)
library(tidyverse)
library(plyr)
library(dplyr)
library(data.table)
library(readr)

setwd("/project/harrislab/naive_vs_infected_whole_brain/R") #run in console

#######################################################################################################

# SPECIFY FILES ON RIVANNA
sample_files <- c( "N1_salmon_quant", 
                   "N2_salmon_quant",
                   "N3_salmon_quant", 
                   "N4_salmon_quant",
                   "I1_salmon_quant",  
                   "I2_salmon_quant", 
                   "I3_salmon_quant", 
                   "I4_salmon_quant")

names(sample_files) <- c("Naive_1", "Naive_2", "Naive_3", "Naive_4", "Infected_1", "Infected_2", "Infected_3", "Infected_4")


sample_files

#######################################################################################################

# SET REFERENCE GENOME

listEnsembl()

ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl", host = "https://uswest.ensembl.org")

#attributes = listAttributes(ensembl)

ensembl_mm <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id", "external_gene_name"), 
                    mart = ensembl) #This is what we need for tximport

colnames(ensembl_mm) = c("TXNAME", "GENEID", "GENENAME")

tx2gene <- ensembl_mm[,1:2] #include only Txname and gene ID


#######################################################################################################

# IMPORT FILES WITH TXIMPORT

setwd("/project/harrislab/naive_vs_infected_whole_brain/salmon/aligned_reads/quant_files")

## Transcript-level 
txi_transcript_level <- tximport(files = sample_files, 
                                 type = "salmon", 
                                 tx2gene = tx2gene,
                                 txOut = TRUE, #set to TRUE to generate the transcript-level output
                                 ignoreTxVersion = TRUE) #set to false if you care about transcript versions/alternative splicing
## Gene-level
txi_gene_level <- tximport(files = sample_files, 
                           type = "salmon", 
                           tx2gene = tx2gene,
                           txOut = FALSE, 
                           ignoreTxVersion = TRUE)


## Make rownames of your table align with the colnames of txi$counts
sampleTable <- data.frame(condition = factor(c("Naive", "Naive", "Naive", "Naive", "Infected", "Infected", "Infected", "Infected")))

rownames(sampleTable) <- colnames(txi_transcript_level$counts) #this is your colData for DESeq2, raw counts
rownames(sampleTable) <- colnames(txi_gene_level$counts)

colnames(txi_transcript_level$counts)
colnames(txi_gene_level$counts)

#txi_transcript_level_df <- as.data.frame(txi_transcript_level) #if you want to manipulate the txi object as a df
#txi_gene_level_df <- as.data.frame(txi_gene_level)

#######################################################################################################

# DIFFERENTIAL EXPRESSION WITH DESeq2

dds_g_raw_counts <- DESeqDataSetFromTximport(txi_gene_level, sampleTable, ~condition)
dds_g_raw_counts <- DESeq(dds_g_raw_counts) ##Performs differential expression analysis based on negative binomial (gamma-poisson) distribution, with default arguments

dds_t_raw_counts <- DESeqDataSetFromTximport(txi_transcript_level, sampleTable, ~condition)
dds_t_raw_counts <- DESeq(dds_t_raw_counts)

## Generate results
res_g <- results(dds_g_raw_counts, contrast = c("condition","Infected", "Naive")) ## define numerator and denominator. set the group you want on the left of a volcano plot first. 
res_t <- results(dds_t_raw_counts, contrast = c("condition","Infected", "Naive"))

summary(res_g)
summary(res_t)

## Save as dataframe
res_g <- as.data.frame(res_g)
res_t <- as.data.frame(res_t)

# For transcript-level analysis, remove transcripts with baseMean = 0 to reduce file size
res_t_filtered <- subset(res_t, baseMean != 0) 

#######################################################################################################

# ANNOTATE RESULTS WITH BIOMART

## Remove row names
res_g_names <- tibble::rownames_to_column(res_g, "ensembl_gene_id")
res_t_names <- tibble::rownames_to_column(res_t_filtered, "ensembl_transcript_id") #note this is on filtered results

## Remove transcript id versions (keep for alternative splicing)
res_t_names$ensembl_transcript_id <- gsub(res_t_names$ensembl_transcript_id, pattern="\\.[[:digit:]]+$", replacement="")

## Fetch matrix with gene names and ENSEMBL gene IDs
### Gene-level
ensembl <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
ensembl_mm_gene <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"), 
                         mart = ensembl)

### Transcript-level
ensembl_mm_transcript <- getBM(attributes = c("external_gene_name", "ensembl_transcript_id"), 
                               mart = ensembl)

## Annotate results with left_join()
res_g_indexed <- left_join(res_g_names, ensembl_mm_gene, by = "ensembl_gene_id")
res_g_indexed <- res_g_indexed[,c(8, 1, 2, 3, 4, 5, 6, 7)] #reordered to preference

res_t_indexed <- left_join(res_t_names, ensembl_mm_transcript, by = "ensembl_transcript_id")
res_t_indexed <- res_t_indexed[,c(8, 1, 2, 3, 4, 5, 6, 7)] #reordered to preference

#######################################################################################################

# EXTRACT NORMALIZED DATA FROM DESEQ2 RESULTS

## Gene-level
normalized_counts_g <- counts(dds_g_raw_counts, normalized = TRUE)
normalized_counts_g <- as.data.frame(normalized_counts_g)
normalized_counts_g <- tibble::rownames_to_column(normalized_counts_g, "ensembl_gene_id")
normalized_counts_g <- left_join(normalized_counts_g, ensembl_mm_gene, by = "ensembl_gene_id")
normalized_counts_g <- normalized_counts_g[,c(10, 1, 2, 3, 4, 5, 6, 7, 8, 9)]


## Transcript-level
normalized_counts_t <- counts(dds_t_raw_counts, normalized = TRUE)
normalized_counts_t <- as.data.frame(normalized_counts_t)
normalized_counts_t <- tibble::rownames_to_column(normalized_counts_t, "ensembl_transcript_id")

normalized_counts_t$ensembl_transcript_id <- gsub(normalized_counts_t$ensembl_transcript_id, pattern="\\.[[:digit:]]+$", replacement="") #remove transcript versions

normalized_counts_t <- left_join(normalized_counts_t, ensembl_mm_transcript, by = "ensembl_transcript_id")
normalized_counts_t <- normalized_counts_t[,c(10, 1, 2, 3, 4, 5, 6, 7, 8, 9)]

#######################################################################################################

# APPLY VARIANCE STABILIZING TRANSFORMATION (VST) 

## Gene-level
vsd_g <- vst(dds_g_raw_counts) 
write.csv(assay(vsd_g), file="gene_level_vst_counts.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
vsd_g_df <- read.csv('gene_level_vst_counts.csv', sep = ",")

names(vsd_g_df)[names(vsd_g_df) == 'X'] <- 'ensembl_gene_id' 

vsd_g_df <- left_join(vsd_g_df, ensembl_mm_gene, by = "ensembl_gene_id")
vsd_g_df <- vsd_g_df[,c(10, 1:9)]


## Transcript-level
vsd_t <- vst(dds_t_raw_counts) 
write.csv(assay(vsd_t), file="transcript_level_vst_counts.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
vsd_t_df <- read.csv('transcript_level_vst_counts.csv', sep = ",")

names(vsd_t_df)[names(vsd_t_df) == 'X'] <- 'ensembl_transcript_id' 

vsd_t_df$ensembl_transcript_id <- gsub(vsd_t_df$ensembl_transcript_id, pattern="\\.[[:digit:]]+$", replacement="")

vsd_t_df <- left_join(vsd_t_df, ensembl_mm_transcript, by = "ensembl_transcript_id")
vsd_t_df <- vsd_t_df[,c(10, 1:9)]

#######################################################################################################

# APPLY RLOG TRANSFORMATION

## Gene-level
rld_g <- rlog(dds_g_raw_counts) 
write.csv(assay(rld_g), file="gene_level_rlog_counts.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
rld_g_df <- read.csv('gene_level_rlog_counts.csv', sep = ",")

names(rld_g_df)[names(rld_g_df) == 'X'] <- 'ensembl_gene_id' 

rld_g_df <- left_join(rld_g_df, ensembl_mm_gene, by = "ensembl_gene_id")
rld_g_df <- rld_g_df[,c(10, 1:9)]


## Transcript-level
rld_t <- rlog(dds_t_raw_counts) 
write.csv(assay(rld_t), file="transcript_level_rlog_counts.csv", quote = FALSE, row.names = TRUE, col.names = TRUE)
rld_t_df <- read.csv('transcript_level_rlog_counts.csv', sep = ",")

names(rld_t_df)[names(rld_t_df) == 'X'] <- 'ensembl_transcript_id' 

rld_t_df$ensembl_transcript_id <- gsub(rld_t_df$ensembl_transcript_id, pattern="\\.[[:digit:]]+$", replacement="")

rld_t_df <- left_join(rld_t_df, ensembl_mm_transcript, by = "ensembl_transcript_id")
rld_t_df <- rld_t_df[,c(10, 1:9)]

#######################################################################################################

# PRINCIPAL COMPONENT ANALYSIS

## Basic visualization
mypca_rld<- plotPCA(rld_g, intgroup = "condition")
mypca_rld

## Customized visualization
pcaData <- plotPCA(rld_g, intgroup= "condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pca2 <- ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=5, shape = 1, stroke = 2) +
  scale_color_manual(values = c("Naive" = "black", "Infected" = "red")) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_minimal() +
  ylim(-20, 20) +
  xlim(-40,40) 
pca2

#######################################################################################################

# CALCULATE FPKM MATRIX

## Gene-level
gene_level_fpkm <- as.data.frame(fpkm(dds_g_raw_counts, robust = TRUE))

### Add columns for group averages
gene_level_fpkm$Naive_FPKM_Mean <- rowMeans(gene_level_fpkm[,1:4])
gene_level_fpkm$Infected_FPKM_mean <- rowMeans(gene_level_fpkm[,5:8])

### Add columns for group standard deviations
gene_level_fpkm$Naive_FPKM_SD <- apply(gene_level_fpkm[,1:4], 1, sd, na.rm = FALSE)
gene_level_fpkm$Infected_FPKM_SD <- apply(gene_level_fpkm[,5:8], 1, sd, na.rm = FALSE)

### Rownames to columns for left join
gene_level_fpkm <- tibble::rownames_to_column(gene_level_fpkm, "ensembl_gene_id")

### Add FPKM statistics to results table
gene_level_fpkm_stats <- gene_level_fpkm[,c(1, 10, 11, 12, 13)] #include only averages and SD columns
res_g_indexed_with_fpkm <- left_join(res_g_indexed, gene_level_fpkm_stats, by = "ensembl_gene_id")
head(res_g_indexed_with_fpkm)


######################################

## Transcript-level
transcript_level_fpkm <- as.data.frame(fpkm(dds_t_raw_counts, robust = TRUE))

### Add columns for group averages
transcript_level_fpkm$Naive_FPKM_Mean <- rowMeans(transcript_level_fpkm[,1:4])
transcript_level_fpkm$Infected_FPKM_Mean <- rowMeans(transcript_level_fpkm[,5:8])

### Add columns for group standard deviations
transcript_level_fpkm$Naive_FPKM_SD <- apply(transcript_level_fpkm[,1:4], 1, sd, na.rm = FALSE)
transcript_level_fpkm$Infected_FPKM_SD <- apply(transcript_level_fpkm[,5:8], 1, sd, na.rm = FALSE)

### Rownames to columns for left join
transcript_level_fpkm <- tibble::rownames_to_column(transcript_level_fpkm, "ensembl_transcript_id")
transcript_level_fpkm$ensembl_transcript_id <- gsub(transcript_level_fpkm$ensembl_transcript_id, pattern="\\.[[:digit:]]+$", replacement="")

### Add FPKM statistics to results table
transcript_level_fpkm_stats <- transcript_level_fpkm[,c(1, 10, 11, 12, 13)] #include only averages and SD columns
res_t_indexed_with_fpkm <- left_join(res_t_indexed, transcript_level_fpkm_stats, by = "ensembl_transcript_id")
head(res_t_indexed_with_fpkm)

######################################

## Annotate FPKM raw counts matrixes
gene_level_fpkm <- left_join(gene_level_fpkm, ensembl_mm_gene, by = "ensembl_gene_id") #add gene names
gene_level_fpkm <- gene_level_fpkm[,c(14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)] #reorder columns

transcript_level_fpkm <- left_join(transcript_level_fpkm, ensembl_mm_transcript, by = "ensembl_transcript_id")
transcript_level_fpkm <- transcript_level_fpkm[,c(14, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)]

#######################################################################################################

# GENE ONTOLOGY ANALYSIS

selection <- function(allScore){ return(allScore < 0.05)} # function that returns TRUE/FALSE for p-values<0.05

## Subset results table to generate a named vector of genes and p-values (input for topGO)
## Genes enriched in Infected relative to Naive
up_genes <- res_g_indexed[res_g_indexed$log2FoldChange > 0, ]
up_genes <- na.omit(up_genes)

## Genes enriched in Naive relative to Infected
down_genes <- res_g_indexed[res_g_indexed$log2FoldChange < 0, ]
down_genes <- na.omit(down_genes)

## Set topGO with appropriate gene universe

### GO Up (enriched in infected group, positive log2FC)
topGO_pvalues_up <- as.numeric(up_genes[, 8])
head(topGO_pvalues_up)

names(topGO_pvalues_up) <- up_genes[, 2]
head(topGO_pvalues_up)

genes_up <- topGO_pvalues_up
genes_up <- na.omit(genes_up)


## GO Down (enriched in naive group, negative log2FC)
topGO_pvalues_down <- as.numeric(down_genes[, 8])
head(topGO_pvalues_down)

names(topGO_pvalues_down) <- down_genes[, 2]
head(topGO_pvalues_down)

genes_down <- topGO_pvalues_down
genes_down <- na.omit(genes_down)

######################################

## Run GO for Biological Processes
allGO2genes <- annFUN.org(whichOnto="BP", feasibleGenes=NULL, mapping="org.Mm.eg.db", ID="ENSEMBL")

### Up genes
GOdata_up <- new("topGOdata",
                 ontology="BP",
                 allGenes= genes_up,
                 annot=annFUN.GO2genes,
                 GO2genes=allGO2genes,
                 geneSel=selection,
                 nodeSize=10) #this argument prunes the GO hierarchy from the terms that have less than 10 annotated genes

### Down genes
GOdata_down <- new("topGOdata",
                   ontology="BP",
                   allGenes= genes_down,
                   annot=annFUN.GO2genes,
                   GO2genes=allGO2genes,
                   geneSel=selection,
                   nodeSize=10)



## Generate GO Table
#runTest is used to apply the specified test statistic and method to the data. 
results.ks_up <- runTest(GOdata_up, algorithm="classic", statistic="ks") #number of nontrivial nodes: 4862
results.ks_down <- runTest(GOdata_down, algorithm="classic", statistic="ks") #number of nontrivial nodes: 4456

goEnrichment_up <- GenTable(GOdata_up, KS=results.ks_up, orderBy="KS", topNodes=4862)
goEnrichment_down <- GenTable(GOdata_down, KS=results.ks_down, orderBy="KS", topNodes=4456)

#######################################################################################################

# EXPORT RESULTS

## Differential expression results
write.table(res_g_indexed_with_fpkm, '/project/harrislab/naive_vs_infected_whole_brain/R/gene_level_results.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(res_t_indexed_with_fpkm, '/project/harrislab/naive_vs_infected_whole_brain/R/transcript_level_results.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

## FPKM matrices
write.table(gene_level_fpkm, '/project/harrislab/naive_vs_infected_whole_brain/R/fpkm_gene_matrix.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(transcript_level_fpkm, '/project/harrislab/naive_vs_infected_whole_brain/R/fpkm_transcript_matrix.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

## Gene ontology results
write.table(goEnrichment_up, '/project/harrislab/naive_vs_infected_whole_brain/R/GO_up.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(goEnrichment_down, '/project/harrislab/naive_vs_infected_whole_brain/R/GO_down.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

## Normalized counts
write.table(normalized_counts_g, '/project/harrislab/naive_vs_infected_whole_brain/R/gene_level_normalized_counts.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(normalized_counts_t, '/project/harrislab/naive_vs_infected_whole_brain/R/transcript_level_normalized_counts.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

## Variance stabilized transformation matrices
write.table(vsd_g_df, '/project/harrislab/naive_vs_infected_whole_brain/R/gene_level_vst_counts.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(vsd_t_df, '/project/harrislab/naive_vs_infected_whole_brain/R/transcript_level_vst_counts.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

## rlog-transformed matrices
write.table(rld_g_df, '/project/harrislab/naive_vs_infected_whole_brain/R/gene_level_rlog_counts.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
write.table(rld_t_df, '/project/harrislab/naive_vs_infected_whole_brain/R/transcript_level_rlog_counts.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

#######################################################################################################