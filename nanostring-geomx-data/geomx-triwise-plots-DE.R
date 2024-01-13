# TRIWISE PLOTS FOR NANOSTRING DATA

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

# BULK OLFACTORY EPITHELIUM

## Subset on 'target_data'
Eoi_replicates_OE <- target_data[, phenoData(target_data)$region==("Bulk olfactory epithelium")]
dim(Eoi_replicates_OE)

## Extract and log2-transform count data
Eoi_OE <- limma::avearrays(Eoi_replicates_OE, phenoData(Eoi_replicates_OE)$class)

Eoi_OE <- log2(Eoi_OE)
Eoi_OE <- Eoi_OE[,c(2,1,3)] # Re-order columns for plot order

dim(Eoi_OE)
colnames(Eoi_OE)

## Perform differential gene expression
design = as(phenoData(Eoi_replicates_OE), "data.frame")
design$class = factor(as.character(design$class))

design <- model.matrix(~0+class, design)
fit <- lmFit(Eoi_replicates_OE, design)
fit = contrasts.fit(fit, matrix(c(1, -1, 0, 0, 1, -1, 1, 0, -1), ncol=3))
fit = eBayes(fit)


## Create a table of top differentially expressed genes
top = topTable(fit, adjust.method = "BH", p.value=0.05, number=Inf, sort.by = "F")

## set list of diffrentially expressed genes to rownames of generated table
Gdiffexp_OE = rownames(top)
summary(Gdiffexp_OE)

## All significant genes with Log2FC > 1 and adjusted p-value of 0.05
Gdiffexp_OE = topTable(fit, adjust.method = "BH", p.value = 0.05, lfc = 1, number = Inf, sort.by = "none")
Gdiffexp_OE = row.names(Gdiffexp_OE)

## Pinned genes 
#Gpin =  topTable(fit, adjust.method = "BH", p.value=0.05, lfc = 1, number=30, sort.by = "F")
#Gpin = rownames(Gpin)
#Gpin = c(Gpin, "Ifnl2", "Ugt2a1")
#summary(Gpin)

## Set genes of interest for each condition
mock_enriched_OE = c("Ugt2a1", "S100a5", "Glcci1", "Rnasek", "Cyp2f2", "Abca13", "Pcp4l1")
virus_enriched_OE = c("Ifnl2", "Egr1", "Junb", "Cxcl10", "Irf1", "Ctss", "Trafd1", "Ifi44", "Ogfr", "Cd274", "Tap1", "Bst2", "Ifit1", "Ifit3", "Isg15", "Parp14", "Usp18", "Serpina3n", "B2m", "Lgals3bp", "H2-K1", "Irf9", "Irf7")

Gpin_OE = c(mock_enriched_OE, virus_enriched_OE)

## Top 20 differentially expressed genes are in red
top20_OE =  topTable(fit, adjust.method = "BH", p.value=0.05, lfc = 1, number=20, sort.by = "F")
top20_OE = rownames(top20_OE)
summary(top20_OE)


## Calculate barycentric coordinates with x and y in separate columns using triwise function
barycoords_OE = transformBarycentric(Eoi_OE)
str(barycoords_OE)

## Dotplot
plotDotplot(barycoords_OE)

## Interactive Dotplot
#install.packages("htmlwidgets")
library(htmlwidgets)

dotplot_OE = interactiveDotplot(Eoi_OE,
                                Gpin = c(Gpin_OE, "S100a4", "Pcp4l1", "Zbp1"),
                                rmax = 5, 
                                Goi = top20_OE
)
dotplot_OE

#saveWidget(dotplot_OE, file = "OE_triwise.html")

#######################################################################################################

# OLFACTORY SENSORY NEURONS

## Subset on OSNs
Eoi_replicates_OSN <- target_data[, phenoData(target_data)$region==("Olfactory sensory neurons")]
dim(Eoi_replicates_OSN)


## Extract and log2-transform count data
Eoi_OSN <- limma::avearrays(Eoi_replicates_OSN, phenoData(Eoi_replicates_OSN)$class)

Eoi_OSN <- log2(Eoi_OSN)
Eoi_OSN <- Eoi_OSN[,c(3,1,2)] #re-order columns for triwise plot order
dim(Eoi_OSN)
colnames(Eoi_OSN)

## Perform Differential Expression Analysis
design = as(phenoData(Eoi_replicates_OSN), "data.frame")
design$class = factor(as.character(design$class))

design <- model.matrix(~0+class, design)
fit <- lmFit(Eoi_replicates_OSN, design)
fit = contrasts.fit(fit, matrix(c(1, -1, 0, 0, 1, -1, 1, 0, -1), ncol=3))

fit = eBayes(fit)


## Select top significant genes
top = topTable(fit, adjust.method = "BH", p.value=0.05, number=Inf, sort.by = "F")

Gdiffexp_OSN = topTable(fit, adjust.method = "BH", p.value = 0.05, lfc = 1, number = Inf, sort.by = "F")
Gdiffexp_OSN = row.names(Gdiffexp_OSN)

### Define genes of interest for eachc group for Gpin argument
mock_enriched = c("Rap1gds1", "S100a5", "Pcp4l1", "Olfr1507", "Rgs19", "Olfr70", "Olfr272")
virus_enriched = c("Trafd1", "Cd274", "Bst2", "Dtx3l", "Oas1a", "Ogfr", "Parp14", "Tapbp", "Rnf213", "Isg15", "Lgals3bp", "B2m", "Ifitm3", "Usp18", "Ifit1", "Ddx60", "Shisa5", "Oasl2", "Irf7", "Irf9", "Rtp4")
Gpin_OSN = c(mock_enriched, virus_enriched)

### Top 20 differentially expressed in red
Goi_OSN =  topTable(fit, adjust.method = "BH", p.value=0.05, lfc = 1, number=20, sort.by = "F")
Goi_OSN = rownames(Goi_OSN)
summary(Goi_OSN)


## Interactive html widget
dotplot = interactiveDotplot(Eoi_OSN, 
                             Gpin = c(Gpin_OSN, "S100a5", "Pcp4l1", "Zbp1"),
                             rmax = 5, 
                             Goi = Goi_OSN
)
dotplot

#saveWidget(dotplot, file = "OSN_triwise.html")

#######################################################################################################
