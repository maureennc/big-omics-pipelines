# Set working directory

dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/mg-sting-ko/results/de'
setwd(dir)
getwd()

library(readr)
library(EnhancedVolcano)

##############################################################################

# IMPORT DATA

library(readr)

cell_type <- read_csv('mg-sting-wt-ko-DE-cell_type-median-norm.csv')
#cell_type$log2fc <- -cell_type$log2fc

#cell_type$log2FoldChange <- as.numeric(cell_type$log2fc)
#cell_type$pval <- as.numeric(cell_type$pval)
#cell_type$qval <- as.numeric(cell_type$qval)
#cell_type$qval <- as.numeric(as.character(cell_type$qval))

## Create data subsets
microglia <- cell_type[cell_type$cell_type == 'Microglia', ]
astrocyte <- cell_type[cell_type$cell_type == 'Astrocyte', ]
excitatory <- cell_type[cell_type$cell_type == 'Excitatory neuron', ]
inhibitory <- cell_type[cell_type$cell_type == 'Inhibitory neuron', ]
oligodendrocyte <- cell_type[cell_type$cell_type == 'Oligodendrocyte', ]
opc <- cell_type[cell_type$cell_type == 'OPC', ]
endothelial <- cell_type[cell_type$cell_type == 'Endothelial cell', ]
pericyte <- cell_type[cell_type$cell_type == 'Pericyte', ]


microglia_filtered <- microglia[microglia$qval < 0.05 & abs(microglia$log2fc) >= 0.5, ]
astrocyte_filtered <- astrocyte[astrocyte$qval < 0.05 & abs(astrocyte$log2fc) >= 0.5, ]
excitatory_filtered <- excitatory[excitatory$qval < 0.05 & abs(excitatory$log2fc) >= 0.5, ]
inhibitory_filtered <- inhibitory[inhibitory$qval < 0.05 & abs(inhibitory$log2fc) >= 0.5, ]
oligodendrocyte_filtered <- oligodendrocyte[oligodendrocyte$qval < 0.05 & abs(oligodendrocyte$log2fc) >= 0.5, ]
opc_filtered <- opc[opc$qval < 0.05 & abs(opc$log2fc) >= 0.5, ]
endothelial_filtered <- endothelial[endothelial$qval < 0.05 & abs(endothelial$log2fc) >= 0.5, ]
pericyte_filtered <- pericyte[pericyte$qval < 0.05 & abs(pericyte$log2fc) >= 0.5, ]

##############################################################################

# VOLCANO PLOT

## Microglia

microglia_genes = c('Tmem173', 'Gm4951', 'Tyk2', 'Kcnt2', 'mt-Nd3', 'Ell2', 'Bcr',
                    'Ttr', 'Hs6st3', 'Zmat4', 'Dkk3', 'Irf9', 'Npas3', 'Fam135b', 'Plp1', 'Ttr')

p1 <- EnhancedVolcano(microglia,
                      lab = microglia$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Microglia',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = microglia_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      #xlim = c(-2, 2),
                      ylim = c(0, 12),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
)
p1

##############################################################################

## Excitatory neuron

excitatory_genes = c('St18', 'Mal', 'Cldn11', 'Ror1', ' Stard13', 'Cd9', 'Adarb2',
                    'Apod', 'Ttr', 'Arhgap6')


p2 <- EnhancedVolcano(excitatory,
                      lab = excitatory$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Excitatory neurons',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = excitatory_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      #xlim = c(-2, 2),
                      #ylim = c(0, 10),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
)

p2

##############################################################################

## Inhibitory neuron

inhibitory_genes = c('Mal', 'Casz1', 'Sox5os4', 'Scn5a', 'Vip', 'Drd3', 'Mdga1', 'Ttr')

p3 <- EnhancedVolcano(inhibitory,
                      lab = inhibitory$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Inhibitory neurons',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = inhibitory_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      xlim = c(-2, 2),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
                      
)

p3

##############################################################################

astrocyte_genes = c('St18', 'Spock3', 'C1ql3', 'Zfp804a', 'Ano3', 'Fmn1', 'Mal', 'Ugt8a', 'C1qa', 'Sel1l3', 'Rimbp2', 'Zfp536', 'Ttr',
                    'Fkbp5', 'Gm20754', 'Cx3cl1', 'Socs2')

p4 <- EnhancedVolcano(astrocyte,
                      lab = astrocyte$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Astrocytes',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = astrocyte_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      #xlim = c(-2, 2),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
                      
)

p4


##############################################################################

## Oligodendrocytes

oligo_genes = c('Rgs20', 'Gli2', 'Ryr3', 'Gm3764', 'Ndst3', 'Ak5', 'Sel1l3', 'Vgf', 'Hif3a', 'Ctsd', 'Nell1', 'Ttr', 'Alk',
                'Grn', 'Alk', 'Hs6st3', 'Lingo1', 'Sv2b', 'Slit3')

p5 <- EnhancedVolcano(oligodendrocyte,
                      lab = oligodendrocyte$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Oligodendrocytes',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = oligo_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      #xlim = c(-2, 2),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
                      
)

p5

##############################################################################

## OPC

opc_genes = c('Cst3', 'Wwtr1', 'Ccdc114', 'A230006K03Rik', 'Ddx60', 'Smpd3', 'Hcn1', 'Airn')

p6 <- EnhancedVolcano(opc,
                      lab = opc$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'OPC',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = opc_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      xlim = c(-2, 2),
                      #ylim = c(0, 17),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
                      
)

p6

##############################################################################

## Endothelial cells


endothelial_genes = c('Fmo2', 'Pde4b', 'C1qa', 'Kcnc2', 'Nrn1', 'Hs6st3', 'Ly6c1', 'Filip1l',
                      'Epha6', 'Gm49678', 'Pde10a', 'Ccnd3', 'Ttr', 'Dcc', 'Myof')

p7 <- EnhancedVolcano(endothelial,
                      lab = endothelial$gene,
                      x = 'log2fc', 
                      y = 'pval', 
                      pCutoffCol = 'qval',
                      title = 'Endothelial cells',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = endothelial_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      xlim = c(-2, 2),
                      #ylim = c(0, 6),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
                      
)

p7

##############################################################################

## Pericytes

pericyte_genes = c('Ttr', 'Ctsd', 'Hcn1', 'Tiam1', 'Mal', 'C1qc',
                   'Tenm1', 'Adgrl4', 'Gabrg3', 'Egfem1',
                   'Strip2', 'Cldn11', 'Smoc2', 'Slc8a1', 'Plp1')

p8 <- EnhancedVolcano(pericyte,
                      lab = pericyte$gene,
                      x = 'log2fc', 
                      y = 'qval', 
                      pCutoffCol = 'qval',
                      title = 'Pericytes',
                      pCutoff = 0.05, 
                      FCcutoff = 0.5,
                      pointSize = 3.0,
                      labSize = 7,
                      selectLab = pericyte_genes,
                      #boxedLabels = TRUE,
                      drawConnectors = TRUE,
                      widthConnectors = 0.75,
                      xlim = c(-2, 2),
                      ylim = c(0, 8),
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      #col=c('grey', 'grey', 'grey', 'red'),
                      legendLabSize = 10,
                      max.overlaps = Inf
                      
)

p8


##############################################################################

# Save the plots
save_dir = '/Users/maureen/Documents/projects/lukens-lab/jess/Thanos_AD-Sting_snRNAseq/analysis-4/mg-sting-ko/results/de/volcano'
setwd(save_dir)

plot_list <- list(p1, p2, p3, p4, p5, p6, p7, p8)

plot_names <- c("Microglia", "Excitatory_neuron", "Inhibitory_neuron",
                "Astrocyte", "Oligodendrocyte", "OPC", "Endothelial_cell", "Pericyte")

## Save as pdf
for (i in seq_along(plot_list)) {
  ggsave(filename = paste0(plot_names[i], ".pdf"),
         plot = plot_list[[i]],
         width = 6, height = 6, dpi = 300,
         path = getwd())
}

## Save as png
for (i in seq_along(plot_list)) {
  ggsave(filename = paste0(plot_names[i], ".png"),
         plot = plot_list[[i]],
         width = 6, height = 6, dpi = 300,
         path = getwd())
}


##############################################################################
