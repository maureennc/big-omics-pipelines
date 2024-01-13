# TRIWISE PLOTS FOR SARS-CoV-2 GENE ONTOLOGY (Bulk Olfactory Epithelium)

## Documentation: https://zouter.github.io/triwise/rd.html#plotdotplot
## Vignette: https://zouter.github.io/triwise/vignette.html

## Run after 'analyze-geomx-data.R'

#######################################################################################################

# PACKAGES

library(triwise)
library(Biobase)
library(limma)
library(ggplot2)
library(org.Mm.eg.db)
library(GO.db)
library(dplyr)

#######################################################################################################

# Load data from NanoString GeoMx Set
Eoi_replicates <- target_data[, phenoData(target_data)$region==("Bulk olfactory epithelium")]
dim(Eoi_replicates)

# Extract count_data
Eoi <- limma::avearrays(Eoi_replicates, phenoData(Eoi_replicates)$class)

#log2 Transform data
Eoi <- log2(Eoi)
Eoi <- Eoi[,c(2,1,3)] #re-order columns for triwise plot order
dim(Eoi)
colnames(Eoi)

# Compute barycentric coordinates
barycoords = transformBarycentric(Eoi)
str(barycoords)

# Perform differential expression
design = as(phenoData(Eoi_replicates_OE), "data.frame")
design$class = factor(as.character(design$class))

design <- model.matrix(~0+class, design)
fit <- lmFit(Eoi_replicates_OE, design)
fit = contrasts.fit(fit, matrix(c(1, -1, 0, 0, 1, -1, 1, 0, -1), ncol=3))
fit = eBayes(fit)

# Top significant genes
top = topTable(fit, adjust.method = "BH", p.value=0.05, number=Inf, sort.by = "F")

Gdiffexp = rownames(top)
summary(Gdiffexp)

# Basic plotting
## Dotplot
plotDotplot(barycoords, Gdiffexp)

## Rose plot
plotRoseplot(barycoords, Gdiffexp, relative = F)

# Extract gene sets
gsets = AnnotationDbi::as.list(org.Mm.egGO2ALLEGS)

# Convert Entrez IDs in Gene sets to gene symbols
gsets_symbols <- sapply(gsets, function(gset) {
  gene_symbols <- mapIds(org.Mm.eg.db, 
                         keys=gset, 
                         column="SYMBOL", 
                         keytype="ENTREZID", 
                         multiVals="first")
  return(gene_symbols)
})

# Run the intersection between gsets and Eoi
gsets <- sapply(gsets_symbols, function(gset) intersect(rownames(Eoi), unique(as.character(gset))))

# Create a table of gene ontology information
gsetindex = dplyr::bind_rows(lapply(AnnotationDbi::as.list(GOTERM[names(gsets)]), function(goinfo) {
  tibble(name=Term(goinfo), definition=Definition(goinfo), ontology=Ontology(goinfo), gsetid = GOID(goinfo))
}))

# Filter for BP Ontology
gsets = gsets[gsetindex %>% filter(ontology == "BP") %>% .$gsetid]

# Test whether a gene is specifically up-regulated in a particular direction
scores = testUnidirectionality(barycoords, gsets, Gdiffexp, statistic = "rank", mc.cores = 8, nsamples=1e+6)
scores = left_join(scores, gsetindex, by="gsetid") %>% filter(qval < 0.05) %>% arrange(qval, z)

scores = scores[(scores$qval < 0.05) & (scores$z > 0.15), ]
scores$redundancy = estimateRedundancy(scores, gsets, Gdiffexp)

# Update gsets so it contains GO name
for (go_id in names(gsets)) {
  go_name <- gsetindex$name[which(gsetindex$gsetid == go_id)]
  names(gsets)[names(gsets) == go_id] <- go_name
}

################################################################################################################

# GENE ONTOLOGY: PLOT GENE SETS OF INTEREST

go_terms_of_interest = c('response to virus', 
                         'interferon-mediated signaling pathway', 
                         'positive regulation of macromolecule metabolic process',
                         'antigen processing and presentation of peptide antigen',
                         'programmed cell death',
                         'lymphocyte mediated immunity',
                         'carboxylic acid metabolic process',
                         'oxacid metabolic process',
                         'plasma membrane bounded cell projection organization',
                         'cell differentiation')

plots = list()

for (go_term in go_terms_of_interest) {
  # Find the title in gsetindex
  plot_title = gsetindex$name[gsetindex$name == go_term]
  
  
  # Generate the plot
  plots[[go_term]] = plotDotplot(barycoords, Gdiffexp, Goi=gsets[[go_term]], showlabels = FALSE) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::ggtitle(plot_title %>% strwrap(40) %>% paste(collapse="\n")) +
    ggplot2::theme(axis.text=ggplot2::element_text(size=14),
                   plot.title=ggplot2::element_text(size=8, face="bold", hjust = 0.5))
}


plots[['response to virus']]
plots[['interferon-mediated signaling pathway']]
plots[['positive regulation of macromolecule metabolic process']]
plots[['antigen processing and presentation of peptide antigen']]
plots[['lymphocyte mediated immunity']]
plots[['carboxylic acid metabolic process']]
plots[['oxoacid metabolic process']]
plots[['plasma membrane bounded cell projection organization']]
plots[['cell differentiation']]


################################################################################################################

