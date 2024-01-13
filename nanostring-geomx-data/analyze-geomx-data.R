# SARS-CoV2 Nanostring GeoMx Analysis

#######################################################################################################

# PACKAGES

library(NanoStringNCTools)
library(GeomxTools)
library(GeoMxWorkflows)

if(packageVersion("GeomxTools") < "2.1" & 
   packageVersion("GeoMxWorkflows") >= "1.0.1"){
  stop("GeomxTools and Workflow versions do not match. Please use the same version. 
    This workflow is meant to be used with most current version of packages. 
    If you are using an older version of Bioconductor please reinstall GeoMxWorkflows and use vignette(GeoMxWorkflows) instead")
}

library(scales)
library(reshape2)
library(cowplot) 
library(umap)
library(Rtsne)
library(pheatmap)
library(ggrepel) 

#######################################################################################################

# LOAD DATASET

## DCCFiles: expression count data & sequencing quality metadata. This is included in the DCC-20230413 directory
## PKCFiles: This is the Mm_R_NGS_WTA_v1.0.pkc file in the pkc directory
## SampleAnnotationFile: This is the LabWorksheet_V3.xlsx annotation file

DCCFiles <- dir(path = "./DCC-20230413", full.names = TRUE, recursive = TRUE)
PKCFiles <- dir(path = "./pkc", full.names = TRUE, recursive = TRUE)
SampleAnnotationFile <- dir(path = "./annotation", full.names = TRUE, recursive = TRUE)

library(readxl)
annotation <- read_excel("./annotation/LabWorksheet_V3.xlsx", na = "NA")

# Import data using arguments defined as above
data <- readNanoStringGeoMxSet(dccFiles = DCCFiles,
                               pkcFiles = PKCFiles,
                               phenoDataFile = SampleAnnotationFile,
                               phenoDataSheet = "Template", #name of the excel sheet tab in .xls file
                               phenoDataDccColName = "Sample_ID", #name of first column with DSP files
                               protocolDataColNames = c("aoi", "roi"), #aoi = geometric/negative; roi = region of interest
                               experimentDataColNames = c("panel")) #column for NGS Whole Transcritome Atlas


#######################################################################################################

# STUDY DESIGN

## Modules used
library(knitr)
pkcs <- annotation(data)
modules <- gsub(".pkc", "", pkcs)
kable(data.frame(PKCs = pkcs, modules = modules))

## Sample overview
library(dplyr)
library(ggforce)

## Select and display the annotations we want to show
count_mat <- count(pData(data), class, region)
count_mat

#######################################################################################################

# QC AND PRE-PROCESSING 
# Standard pre-processing includes selecting ROI/AOI segments and genes based on quality control (QC) or limit of quantification (LOQ) metrics to normalize the data

## Shift counts to one
data <- shiftCountsOne(data, useDALogic = TRUE)

## SELECT SEGMENT QC
QC_params <-
  list(minSegmentReads = 1000, # Minimum number of reads (1000)
       percentTrimmed = 80,    # Minimum % of reads trimmed (80%)
       percentStitched = 80,   # Minimum % of reads stitched (80%)
       percentAligned = 75,    # Minimum % of reads aligned (75%)
       percentSaturation = 50, # Minimum sequencing saturation (50%)
       minNegativeCount = 0,   # Minimum negative control counts (1)
       maxNTCCount = 1000,     # Maximum counts observed in NTC well (9000)
       minNuclei = 20,         # Minimum # of nuclei estimated (100)
       minArea = 1000)         # Minimum segment area (5000)

data <- setSegmentQCFlags(data, qcCutoffs = QC_params)    

## Collate QC Results to generate data summary
QCResults <- protocolData(data)[["QCFlags"]]
flag_columns <- colnames(QCResults)
QC_Summary <- data.frame(Pass = colSums(!QCResults[, flag_columns]),
                         Warning = colSums(QCResults[, flag_columns]))
QCResults$QCStatus <- apply(QCResults, 1L, function(x) {
  ifelse(sum(x) == 0L, "PASS", "WARNING")
})

QC_Summary["TOTAL FLAGS", ] <-
  c(sum(QCResults[, "QCStatus"] == "PASS"),
    sum(QCResults[, "QCStatus"] == "WARNING"))

QCResults
QC_Summary


## QC visualization function
col_by <- "segment"

QC_histogram <- function(assay_data = NULL,
                         annotation = NULL,
                         fill_by = NULL,
                         thr = NULL,
                         scale_trans = NULL) {
  plt <- ggplot(assay_data,
                aes_string(x = paste0("unlist(`", annotation, "`)"),
                           fill = fill_by)) +
    geom_histogram(bins = 50) +
    geom_vline(xintercept = thr, lty = "dashed", color = "black") +
    theme_bw() + guides(fill = "none") +
    facet_wrap(as.formula(paste("~", fill_by)), nrow = 4) +
    labs(x = annotation, y = "Segments, #", title = annotation)
  if(!is.null(scale_trans)) {
    plt <- plt +
      scale_x_continuous(trans = scale_trans)
  }
  plt
}


## Histograms for each QC Segment
library(ggplot2) #import ggplot2 library

QC_histogram(sData(data), "Trimmed (%)", col_by, 80)
QC_histogram(sData(data), "Stitched (%)", col_by, 80)
QC_histogram(sData(data), "Aligned (%)", col_by, 75)
QC_histogram(sData(data), "Saturated (%)", col_by, 50) + labs(title = "Sequencing Saturation (%)", x = "Sequencing Saturation (%)")
QC_histogram(sData(data), "area", col_by, 1000, scale_trans = "log10")
QC_histogram(sData(data), "nuclei", col_by, 20)

## Calculate the negative geometric means for each module
negativeGeoMeans <- 
  esBy(negativeControlSubset(data), 
       GROUP = "Module", 
       FUN = function(x) { 
         assayDataApply(x, MARGIN = 2, FUN = ngeoMean, elt = "exprs") 
       }) 
protocolData(data)[["NegGeoMean"]] <- negativeGeoMeans

# explicitly copy the Negative geoMeans from sData to pData (phenotypic data layer in data object)
negCols <- paste0("NegGeoMean_", modules)
pData(data)[, negCols] <- sData(data)[["NegGeoMean"]]
for(ann in negCols) {
  plt <- QC_histogram(pData(data), ann, col_by, 2, scale_trans = "log10")
  print(plt)
}


## QC summary
pData(data) <- pData(data)[, !colnames(pData(data)) %in% negCols] # Detach neg_geomean columns ahead of aggregateCounts call

kable(QC_Summary, caption = "QC Summary Table for each Segment") # Create a table tallying how many ROIs passed each QC parameter

### Four samples are flagged due to failure to pass QC


## Filter data NanostringGeoMxSet for segments that passed QC
data <- data[, QCResults$QCStatus == "PASS"]
dim(data) # There are now 138 samples that successfully passed QC

#######################################################################################################

# PROBE QC

## Grubbs outlier test to remove low-performing probe
data <- setBioProbeQCFlags(data, qcCutoffs = list(minProbeRatio = 0.1, percentFailGrubbs = 20), removeLocalOutliers = FALSE)

ProbeQCPassed <- 
  subset(data, 
         fData(data)[["QCFlags"]][,c("LowProbeRatio")] == FALSE &
           fData(data)[["QCFlags"]][,c("GlobalGrubbsOutlier")] == FALSE)
dim(ProbeQCPassed)

data <- ProbeQCPassed 

#######################################################################################################

# CREATE GENE-LEVEL COUNT DATA ('target_data')

## Check how many unique targets the object has
length(unique(featureData(data)[["TargetName"]]))

target_data <- aggregateCounts(data) 
dim(target_data)

exprs(target_data)

#######################################################################################################

# CALCULATE LIMIT OF QUANTIFICATION FOR EACH SEGMENT
cutoff <- 2
minLOQ <- 2

# Calculate LOQ per module tested
LOQ <- data.frame(row.names = colnames(target_data))
for(module in modules) {
  vars <- paste0(c("NegGeoMean_", "NegGeoSD_"),
                 module)
  if(all(vars[1:2] %in% colnames(pData(target_data)))) {
    LOQ[, module] <-
      pmax(minLOQ,
           pData(target_data)[, vars[1]] * 
             pData(target_data)[, vars[2]] ^ cutoff)
  }
}

pData(target_data)$LOQ <- LOQ


#######################################################################################################

# FILTER SEGMENTS AND GENES WITH POOR SIGNAL

LOQ_Mat <- c()
for(module in modules) {
  ind <- fData(target_data)$Module == module
  Mat_i <- t(esApply(target_data[ind, ], MARGIN = 1,
                     FUN = function(x) {
                       x > LOQ[, module]
                     }))
  LOQ_Mat <- rbind(LOQ_Mat, Mat_i)
}
# ensure ordering since this is stored outside of the geomxSet
LOQ_Mat <- LOQ_Mat[fData(target_data)$TargetName, ]


## Save detection rate information to pheno date

pData(target_data)$GenesDetected <- 
  colSums(LOQ_Mat, na.rm = TRUE) 

pData(target_data)$GeneDetectionRate <-
  pData(target_data)$GenesDetected / nrow(target_data)

## Determine detection thresholds: 1%, 5%, 10%, 15%, >15%
pData(target_data)$DetectionThreshold <- 
  cut(pData(target_data)$GeneDetectionRate,
      breaks = c(0, 0.01, 0.05, 0.1, 0.15, 1),
      labels = c("<1%", "1-5%", "5-10%", "10-15%", ">15%"))

## stacked bar plot of different cut points (1%, 5%, 10%, 15%)
ggplot(pData(target_data),
       aes(x = DetectionThreshold)) +
  geom_bar(aes(fill = region)) +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_bw() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = "Gene Detection Rate",
       y = "Segments, #",
       fill = "Segment Type")


## Remove segments with <5% gene detection rate
target_data<-
  target_data[, pData(target_data)$GeneDetectionRate >= .05]

dim(target_data)

#######################################################################################################

# DETERMINE DETECTION RATE FOR GENES

#library(scales)

# Calculate detection rate:
LOQ_Mat <- LOQ_Mat[, colnames(target_data)]
fData(target_data)$DetectedSegments <- rowSums(LOQ_Mat, na.rm = TRUE)
fData(target_data)$DetectionRate <-
  fData(target_data)$DetectedSegments / nrow(pData(target_data))

# Gene of interest detection table
goi <- c("Bst2", "Ugt2a1", "Ugt2a2", "Irf1", "Irf8", "Irf9") # example genes
goi_df <- data.frame(
  Gene = goi,
  Number = fData(target_data)[goi, "DetectedSegments"],
  DetectionRate = percent(fData(target_data)[goi, "DetectionRate"]))

goi_df

#######################################################################################################

# FILTER GENES

## Plot detection rate:
plot_detect <- data.frame(Freq = c(1, 5, 10, 20, 30, 50))
plot_detect$Number <-
  unlist(lapply(c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5),
                function(x) {sum(fData(target_data)$DetectionRate >= x)}))
plot_detect$Rate <- plot_detect$Number / nrow(fData(target_data))
rownames(plot_detect) <- plot_detect$Freq

ggplot(plot_detect, aes(x = as.factor(Freq), y = Rate, fill = Rate)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = formatC(Number, format = "d", big.mark = ",")),
            vjust = 1.6, color = "black", size = 4) +
  scale_fill_gradient2(low = "orange2", mid = "lightblue",
                       high = "dodgerblue3", midpoint = 0.65,
                       limits = c(0,1),
                       labels = scales::percent) +
  theme_bw() +
  scale_y_continuous(labels = scales::percent, limits = c(0,1),
                     expand = expansion(mult = c(0, 0))) +
  labs(x = "% of Segments",
       y = "Genes Detected, % of Panel > LOQ")


## Subset to target genes detected in at least 5% of the samples. Also manually include the negative control probe, for downstream use
negativeProbefData <- subset(fData(target_data), CodeClass == "Negative")
neg_probes <- unique(negativeProbefData$TargetName)
target_data <- 
  target_data[fData(target_data)$DetectionRate >= 0.05 |
                fData(target_data)$TargetName %in% neg_probes, ]
dim(target_data)

## Retain only detected genes of interest
goi <- goi[goi %in% rownames(target_data)]

#######################################################################################################

# Q3 NORMALIZATION

## Graph Q3 Values

#library(reshape2)  # for melt
#library(cowplot)   # for plot_grid

# Graph Q3 value vs negGeoMean of Negatives
ann_of_interest <- "region"
Stat_data <- 
  data.frame(row.names = colnames(exprs(target_data)),
             Segment = colnames(exprs(target_data)),
             Annotation = pData(target_data)[, ann_of_interest],
             Q3 = unlist(apply(exprs(target_data), 2,
                               quantile, 0.75, na.rm = TRUE)),
             NegProbe = exprs(target_data)[neg_probes, ])
Stat_data_m <- melt(Stat_data, measure.vars = c("Q3", "NegProbe"),
                    variable.name = "Statistic", value.name = "Value")

plt1 <- ggplot(Stat_data_m,
               aes(x = Value, fill = Statistic)) +
  geom_histogram(bins = 40) + theme_bw(base_size = 10) +
  scale_x_continuous(trans = "log2") +
  facet_wrap(~Annotation, nrow = 1) + 
  scale_fill_brewer(palette = 3, type = "qual") +
  labs(x = "Counts", y = "Segments, #")

plt2 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3, color = Annotation)) +
  geom_abline(intercept = 0, slope = 1, lty = "dashed", color = "darkgray") +
  geom_point() + guides(color = "none") + theme_minimal(base_size = 10) +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3 Value, Counts")

plt3 <- ggplot(Stat_data,
               aes(x = NegProbe, y = Q3 / NegProbe, color = Annotation)) +
  geom_hline(yintercept = 1, lty = "dashed", color = "darkgray") +
  geom_point() + theme_bw(base_size = 10) +
  scale_x_continuous(trans = "log2") + 
  scale_y_continuous(trans = "log2") +
  theme(aspect.ratio = 1) +
  labs(x = "Negative Probe GeoMean, Counts", y = "Q3/NegProbe Value, Counts")

btm_row <- plot_grid(plt2, plt3, nrow = 1, labels = c("B", ""),
                     rel_widths = c(0.43,0.57))
plot_grid(plt1, btm_row, ncol = 1, labels = c("A", ""))


## Q3 normalization (75th percentile) for WTA/CTA  with or without custom spike-ins
target_data <- normalize(target_data ,
                         norm_method = "quant", 
                         desiredQuantile = .75,
                         toElt = "q_norm")

## Background normalization for WTA/CTA without custom spike-in
target_data <- normalize(target_data ,
                         norm_method = "neg", 
                         fromElt = "exprs",
                         toElt = "neg_norm")



## Boxplots to visualize effects of normalization on the data

### Boxplot for raw counts
boxplot(exprs(target_data)[,1:10],
        col = "#9EDAE5", main = "Raw Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Raw")

## Boxplot for Q3 Normalized counts
boxplot(assayDataElement(target_data[,1:10], elt = "q_norm"),
        col = "#2CA02C", main = "Q3 Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Q3 Normalized")

## Boxplot for Negative Normalized Counts
boxplot(assayDataElement(target_data[,1:10], elt = "neg_norm"),
        col = "#FF7F0E", main = "Neg Norm Counts",
        log = "y", names = 1:10, xlab = "Segment",
        ylab = "Counts, Neg. Normalized")

#######################################################################################################

# DIMENSIONALITY REDUCTION

## UMAP
custom_umap <- umap::umap.defaults
custom_umap$random_state <- 42

umap_out <-
  umap(t(log2(assayDataElement(target_data , elt = "q_norm"))),  
       config = custom_umap)
pData(target_data)[, c("UMAP1", "UMAP2")] <- umap_out$layout[, c(1,2)]


ggplot(pData(target_data),
       aes(x = UMAP1, y = UMAP2, color = class, shape = region)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7))



## t-SNE
set.seed(42)
tsne_out <-
  Rtsne(t(log2(assayDataElement(target_data , elt = "q_norm"))),
        perplexity = ncol(target_data)*.15)
pData(target_data)[, c("tSNE1", "tSNE2")] <- tsne_out$Y[, c(1,2)]


ggplot(pData(target_data),
       aes(x = tSNE1, y = tSNE2, color = class, shape = region)) +
  geom_point(size = 3) +
  theme_classic() +
  theme(panel.border = element_rect(color = "black", fill = NA, size = 0.7))

#######################################################################################################

# CLUSTERING HIGH COEFFICIENT OF VARIATION (CV) GENES

## Log2-transform the data
assayDataElement(object = target_data, elt = "log_q") <-
  assayDataApply(target_data, 2, FUN = log, base = 2, elt = "q_norm")

## Create CV function
calc_CV <- function(x) {sd(x) / mean(x)}
CV_dat <- assayDataApply(target_data,
                         elt = "log_q", MARGIN = 1, calc_CV)
## Show the highest CD genes and their CV values
sort(CV_dat, decreasing = TRUE)[1:5]

## Identify genes in the top 3rd of the CV values
GOI <- names(CV_dat)[CV_dat > quantile(CV_dat, 0.8)]
pheatmap(assayDataElement(target_data[GOI, ], elt = "log_q"),
         scale = "row", 
         show_rownames = FALSE, show_colnames = FALSE,
         border_color = NA,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         breaks = seq(-3, 3, 0.05),
         color = colorRampPalette(c("purple3", "black", "yellow2"))(120),
         annotation_col = 
           pData(target_data)[, c("class","region")])

#######################################################################################################

# DIFFERENTIAL EXPRESSION (BETWEEN-SLIDES; GLOBAL)

## This chunk may take a while to run. It will generate a gene results table for mock vs. SARS-CoV-2 based on region (eptihelium, glomeruli, etc)

## Convert test variables to factors
pData(target_data)$testClass <-
  factor(pData(target_data)$infection, c("SARS-CoV-2", "Mock"))

## Run Linear Mixed Mode; formula follows conventions defined by the lme4 package
results_between <- c()
for(region in c("Bulk olfactory epithelium", "Bulk glomeruli", "Susentacular cells", "Olfactory sensory neurons")) {
  ind <- pData(target_data)$region == region
  mixedOutmc <-
    mixedModelDE(target_data[, ind],
                 elt = "log_q",
                 modelFormula = ~ testClass + (1 | slide), #test
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  # Format results as data.frame
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  # Correctly associate gene name with it's row in the results table
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results_between <- rbind(results_between, r_test)
}

# rownames(results_between) <- c(1:30040)

## Categorize differential expressed genes based on P-value & false discovery rate (FDR) for plotting
results_between$Color <- "NS or FC < 0.5"
results_between$Color[results_between$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results_between$Color[results_between$FDR < 0.05] <- "FDR < 0.05"
results_between$Color[results_between$FDR < 0.001] <- "FDR < 0.001"
results_between$Color[abs(results_between$Estimate) < 0.5] <- "NS or FC < 0.5"
results_between$Color <- factor(results_between$Color,
                                levels = c("NS or FC < 0.5", "P < 0.05",
                                           "FDR < 0.05", "FDR < 0.001"))

results_between$invert_P <- (-log10(results_between$`Pr(>|t|)`)) * sign(results_between$Estimate)

top_g <- c()

for(cond in c("Bulk glomeruli", "Bulk olfactory epithelium", "Susentacular cells", "Olfactory sensory neurons")) {
  ind <- results_between$Subset == cond
  top_g <- c(top_g,
             results_between[ind, 'Gene'][
               order(results_between[ind, 'invert_P'], decreasing = TRUE)[1:20]],
             results_between[ind, 'Gene'][
               order(results_between[ind, 'invert_P'], decreasing = FALSE)[1:20]])
}

top_g <- unique(top_g)

results_between <- results_between[, -1*ncol(results_between)] # remove invert_P from matrix


## Graph results for between-slides comparison
ggplot(results_between,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Mock <- log2(FC) -> SARS-CoV-2",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results_between, Gene %in% top_g & FDR < 0.001),
                  size = 2, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  facet_wrap(~Subset, scales = "free_y")

#######################################################################################################

# DIFFERENTIAL EXPRESSION (BETWEEN-SLIDES, Group A vs. Group B)
## Mock vs. Infected, virus absent

### Subset data
AB <- subset(target_data, select = phenoData(target_data)[["class"]] == c("Mock-infected", "Virus absent"))

### Convert test variables to factors
pData(AB)$testClass <-
  factor(pData(AB)$class, c("Virus absent", "Mock-infected"))

### Run LMM:
results_AB <- c()
for(region in c("Bulk olfactory epithelium", "Bulk glomeruli", "Susentacular cells", "Olfactory sensory neurons")) {
  ind <- pData(AB)$region == region
  mixedOutmc <-
    mixedModelDE(AB[, ind],
                 elt = "log_q",
                 modelFormula = ~ testClass + (1 | slide),
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results_AB <- rbind(results_AB, r_test)
}



rownames(results_AB) <- c(1:33224)


### Categorize Results based on P-value & FDR for plotting
results_AB$Color <- "NS or FC < 0.5"
results_AB$Color[results_AB$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results_AB$Color[results_AB$FDR < 0.05] <- "FDR < 0.05"
results_AB$Color[results_AB$FDR < 0.001] <- "FDR < 0.001"
results_AB$Color[abs(results_AB$Estimate) < 0.5] <- "NS or FC < 0.5"
results_AB$Color <- factor(results_AB$Color,
                           levels = c("NS or FC < 0.5", "P < 0.05",
                                      "FDR < 0.05", "FDR < 0.001"))

### Select top genes
results_AB$invert_P <- (-log10(results_AB$`Pr(>|t|)`)) * sign(results_AB$Estimate)
top_g <- c()
for(cond in c("Mock-infected", "Virus absent")) {
  ind <- results_AB$Subset == cond
  top_g <- c(top_g,
             results_AB[ind, 'Gene'][
               order(results_AB[ind, 'invert_P'], decreasing = TRUE)[1:100]],
             results_AB[ind, 'Gene'][
               order(results_AB[ind, 'invert_P'], decreasing = FALSE)[1:100]])
}
top_g <- unique(top_g)
results_AB <- results_AB[, -1*ncol(results_AB)] # remove invert_P from matrix


### Graph results
ggplot(results_AB,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Mock-infected <- log2(FC) -> Virus absent",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results_AB, Gene %in% top_g & FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 16) +
  theme(legend.position = "bottom") +
  facet_wrap(~Subset, scales = "free_y")
```

#######################################################################################################

# DIFFERENTIAL EXPRESSION (BETWEEN-SLIDES, Group A vs. Group C)
## Mock vs. Infected, virus present (A vs C)

AC <- subset(target_data, select = phenoData(target_data)[["class"]] == c("Mock-infected", "Virus present"))

### Convert test variables to factors
pData(AC)$testClass <-
  factor(pData(AC)$class, c("Virus present", "Mock-infected"))

### Run LMM
results_AC <- c()
for(region in c("Bulk olfactory epithelium", "Susentacular cells", "Olfactory sensory neurons")) {
  ind <- pData(AC)$region == region
  mixedOutmc <-
    mixedModelDE(AC[, ind],
                 elt = "log_q",
                 modelFormula = ~ testClass + (1 | slide),
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests

  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results_AC <- rbind(results_AC, r_test)
}


rownames(results_AC) <- c(1:24918)


## Categorize Results based on P-value & FDR for plotting
results_AC$Color <- "NS or FC < 0.5"
results_AC$Color[results_AC$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results_AC$Color[results_AC$FDR < 0.05] <- "FDR < 0.05"
results_AC$Color[results_AC$FDR < 0.001] <- "FDR < 0.001"
results_AC$Color[abs(results_AC$Estimate) < 0.5] <- "NS or FC < 0.5"
results_AC$Color <- factor(results_AC$Color,
                           levels = c("NS or FC < 0.5", "P < 0.05",
                                      "FDR < 0.05", "FDR < 0.001"))

### Select top genes
results_AC$invert_P <- (-log10(results_AC$`Pr(>|t|)`)) * sign(results_AC$Estimate)
top_g <- c()
for(cond in c("Mock-infected", "Virus present")) {
  ind <- results_AC$Subset == cond
  top_g <- c(top_g,
             results_AC[ind, 'Gene'][
               order(results_AC[ind, 'invert_P'], decreasing = TRUE)[1:100]],
             results_AC[ind, 'Gene'][
               order(results_AC[ind, 'invert_P'], decreasing = FALSE)[1:100]])
}
top_g <- unique(top_g)
results_AC <- results_AC[, -1*ncol(results_AC)] # remove invert_P from matrix

### Graph results
ggplot(results_AC,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Mock-infected <- log2(FC) -> Virus present",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results_AB, Gene %in% top_g & FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  facet_wrap(~Subset, scales = "free_y")

#######################################################################################################

# DIFFERENTIAL EXPRESSION (BETWEEN-SLIDES, Group B vs. Group C)
## Mock vs. Infected, virus absent (B vs. C)

BC <- subset(target_data, select = phenoData(target_data)[["class"]] == c("Virus present", "Virus absent"))

### Convert test variables to factors
pData(BC)$testClass <-
  factor(pData(BC)$class, c("Virus present", "Virus absent"))

### Run LMM
results_BC <- c()
for(region in c("Bulk olfactory epithelium", "Susentacular cells", "Olfactory sensory neurons")) {
  ind <- pData(BC)$region == region
  mixedOutmc <-
    mixedModelDE(BC[, ind],
                 elt = "log_q",
                 modelFormula = ~ testClass + (1 | slide),
                 groupVar = "testClass",
                 nCores = parallel::detectCores(),
                 multiCore = FALSE)
  
  r_test <- do.call(rbind, mixedOutmc["lsmeans", ])
  tests <- rownames(r_test)
  r_test <- as.data.frame(r_test)
  r_test$Contrast <- tests
  
  r_test$Gene <- 
    unlist(lapply(colnames(mixedOutmc),
                  rep, nrow(mixedOutmc["lsmeans", ][[1]])))
  r_test$Subset <- region
  r_test$FDR <- p.adjust(r_test$`Pr(>|t|)`, method = "fdr")
  r_test <- r_test[, c("Gene", "Subset", "Contrast", "Estimate", 
                       "Pr(>|t|)", "FDR")]
  results_BC <- rbind(results_BC, r_test)
}

rownames(results_BC) <- c(1:24918)


### Categorize Results based on P-value & FDR for plotting
results_BC$Color <- "NS or FC < 0.5"
results_BC$Color[results_BC$`Pr(>|t|)` < 0.05] <- "P < 0.05"
results_BC$Color[results_BC$FDR < 0.05] <- "FDR < 0.05"
results_BC$Color[results_BC$FDR < 0.001] <- "FDR < 0.001"
results_BC$Color[abs(results_BC$Estimate) < 0.5] <- "NS or FC < 0.5"
results_BC$Color <- factor(results_BC$Color,
                           levels = c("NS or FC < 0.5", "P < 0.05",
                                      "FDR < 0.05", "FDR < 0.001"))

### Select top genes
results_BC$invert_P <- (-log10(results_BC$`Pr(>|t|)`)) * sign(results_BC$Estimate)
top_g <- c()
for(cond in c("Mock-infected", "Virus absent")) {
  ind <- results_BC$Subset == cond
  top_g <- c(top_g,
             results_BC[ind, 'Gene'][
               order(results_BC[ind, 'invert_P'], decreasing = TRUE)[1:15]],
             results_BC[ind, 'Gene'][
               order(results_BC[ind, 'invert_P'], decreasing = FALSE)[1:15]])
}
top_g <- unique(top_g)
results_BC <- results_BC[, -1*ncol(results_BC)] # remove invert_P from matrix


### Graph results
ggplot(results_BC,
       aes(x = Estimate, y = -log10(`Pr(>|t|)`),
           color = Color, label = Gene)) +
  geom_vline(xintercept = c(0.5, -0.5), lty = "dashed") +
  geom_hline(yintercept = -log10(0.05), lty = "dashed") +
  geom_point() +
  labs(x = "Mock-infected <- log2(FC) -> Virus absent",
       y = "Significance, -log10(P)",
       color = "Significance") +
  scale_color_manual(values = c(`FDR < 0.001` = "dodgerblue",
                                `FDR < 0.05` = "lightblue",
                                `P < 0.05` = "orange2",
                                `NS or FC < 0.5` = "gray"),
                     guide = guide_legend(override.aes = list(size = 4))) +
  scale_y_continuous(expand = expansion(mult = c(0,0.05))) +
  geom_text_repel(data = subset(results_AB, Gene %in% top_g & FDR < 0.05),
                  size = 4, point.padding = 0.15, color = "black",
                  min.segment.length = .1, box.padding = .2, lwd = 2,
                  max.overlaps = 50) +
  theme_bw(base_size = 14) +
  theme(legend.position = "bottom") +
  facet_wrap(~Subset, scales = "free_y")


#######################################################################################################

# EXPORT DIFFERENTIAL EXPRESSION RESULTS

## Prepare results

### A vs. B
AB_OE_bulk <- subset(results_AB, Subset == "Bulk olfactory epithelium")
AB_GL_bulk <- subset(results_AB, Subset == "Bulk glomeruli")
AB_OSN <- subset(results_AB, Subset == "Olfactory sensory neurons")
AB_SUS <- subset(results_AB, Subset == "Susentacular cells")

### A vs. C
AC_OE_bulk <- subset(results_AC, Subset == "Bulk olfactory epithelium")
AC_OSN <- subset(results_AC, Subset == "Olfactory sensory neurons")
AC_SUS <- subset(results_AC, Subset == "Susentacular cells")

### B vs. C
BC_OE_bulk <- subset(results_BC, Subset == "Bulk olfactory epithelium")
BC_OSN <- subset(results_BC, Subset == "Olfactory sensory neurons")
BC_SUS <- subset(results_BC, Subset == "Susentacular cells")


## Write files

### A vs. B
#write.table(AB_OE_bulk, './write/spreadsheets/OE_bulk_Mock_vs_virus_absent.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(AB_GL_bulk, './write/spreadsheets/GL_bulk_Mock_vs_virus_absent.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(AB_OSN, './write/spreadsheets/OSN_Mock_vs_virus_absent.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(AB_SUS, './write/spreadsheets/SUS_Mock_vs_virus_absent.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


### A vs. C
#write.table(AC_OE_bulk, './write/spreadsheets/OE_bulk_Mock_vs_virus_present.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(AC_OSN, './write/spreadsheets/OSN_Mock_vs_virus_present.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(AC_SUS, './write/spreadsheets/SUS_Mock_vs_virus_present.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)


### B vs. C
#write.table(BC_OE_bulk, './write/spreadsheets/OE_bulk_virus_absent_vs_present.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(BC_OSN, './write/spreadsheets/OSN_virus_absent_vs_present.csv.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)
#write.table(BC_SUS, './write/spreadsheets/SUS_virus_absent_vs_present.csv.csv', quote=FALSE, sep=",", row.names=FALSE, col.names=TRUE)

#######################################################################################################
