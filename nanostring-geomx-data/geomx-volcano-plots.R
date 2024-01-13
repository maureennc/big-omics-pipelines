# SARS-CoV-2 GENE EXPRESSION - DATA VISUALIZATION

#######################################################################################################

# LOAD DATA


AB_OE_bulk <- read.csv('./write/spreadsheets/OE_bulk_Mock_vs_virus_absent.csv', header = TRUE, sep = ",")
AB_GL_bulk <- read.csv('./write/spreadsheets/GL_bulk_Mock_vs_virus_absent.csv', header = TRUE, sep = ",")
AB_OSN <- read.csv('./write/spreadsheets/OSN_Mock_vs_virus_absent.csv', header = TRUE, sep = ",")
AB_SUS <- read.csv('./write/spreadsheets/SUS_Mock_vs_virus_absent.csv', header = TRUE, sep = ",")

AC_OE_bulk <- read.csv('./write/spreadsheets/OE_bulk_Mock_vs_virus_present.csv', header = TRUE, sep = ",")
AC_OSN <- read.csv('./write/spreadsheets/OSN_Mock_vs_virus_present.csv', header = TRUE, sep = ",")
AC_SUS <- read.csv('./write/spreadsheets/SUS_Mock_vs_virus_present.csv', header = TRUE, sep = ",")

BC_OE_bulk <- read.csv('./write/spreadsheets/OE_bulk_virus_absent_vs_present.csv', header = TRUE, sep = ",")
BC_OSN <- read.csv('./write/spreadsheets/OSN_virus_absent_vs_present.csv', header = TRUE, sep = ",")
BC_SUS <- read.csv('./write/spreadsheets/SUS_virus_absent_vs_present.csv', header = TRUE, sep = ",")


#######################################################################################################

library(EnhancedVolcano)

#######################################################################################################

# OLFACTORY EPITHELIUM: MOCK VS. VIRUS ABSENT (Group A vs. Group B)

## Set if-else statements for coloring / shapes
keyvals.color <- ifelse(
  AB_OE_bulk$Estimate < 0, '#AFAFB0',
  ifelse(AB_OE_bulk$Estimate > 0, '#FB02FF', 'black'))

names(keyvals.color)[keyvals.color == '#AFAFB0'] <- 'Enriched in Mock-infected'
names(keyvals.color)[keyvals.color == '#FB02FF'] <- 'Enriched in Infected, virus absent'



keyvals.shape <- ifelse(
  AB_OE_bulk$FDR <0.05 & AB_OE_bulk$Estimate < 0, 1,
  ifelse(AB_OE_bulk$FDR < 0.05 & AB_OE_bulk$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


## Plot
p1 <- EnhancedVolcano(AB_OE_bulk,
                      lab = AB_OE_bulk$Gene,
                      #selectLab = c(),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Olfactory Epithelium: Mock vs. Virus Absent',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-2, 3),
                      ylim = c(0, 2),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 20)
p1



#######################################################################################################

# BULK OLFACTORY EPITHELIUM: MOCK VS. VIRUS PRESENT (Group A vs. Group C)

keyvals.color <- ifelse(
  AC_OE_bulk$Estimate < 0, '#AFAFB0',
  ifelse(AC_OE_bulk$Estimate > 0, '#21FF06', 'black'))

names(keyvals.color)[keyvals.color == '#AFAFB0'] <- 'Enriched in Mock-infected'
names(keyvals.color)[keyvals.color == '#21FF06'] <- 'Enriched in Infected, virus present'


keyvals.shape <- ifelse(
  AC_OE_bulk$FDR <0.05 & AC_OE_bulk$Estimate < 0, 1,
  ifelse(AC_OE_bulk$FDR < 0.05 & AC_OE_bulk$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'



p2 <- EnhancedVolcano(AC_OE_bulk,
                      lab = AC_OE_bulk$Gene,
                      selectLab = c("Ugt2a1", "Ugt2a2", "Ifit1", 'Ifit3', "Parp14", "Bst2", 'Isg15', 'Ifitm3', 'Irf9', 'Irf7', 'Cd274', 'H2-T23', 'Tap1', 'Mx1', 'Ifit2', 'Cxcl10', 'Cd74', 'Ccl5', 'Gbp3', 'Stat2', 'Stat1', 'B2m', "Lrp2", "Pclo", "Rnasek", "Abca13", "S100a5"),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Olfactory Epithelium: Mock vs. Virus Present',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-5, 8),
                      ylim = c(0, 8),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 20)
p2


#######################################################################################################

# BULK OLFACTORY EPITHELIUM: VIRUS ABSENT VS. VIRUS PRESENT (Group B vs. Group C)

keyvals.color <- ifelse(
  BC_OE_bulk$Estimate < 0, '#FB02FF',
  ifelse(BC_OE_bulk$Estimate > 0, '#21FF06', 'black'))

names(keyvals.color)[keyvals.color == '#FB02FF'] <- 'Enriched in Infected, virus absent'
names(keyvals.color)[keyvals.color == '#21FF06'] <- 'Enriched in Infected, virus present'


keyvals.shape <- ifelse(
  BC_OE_bulk$FDR <0.05 & BC_OE_bulk$Estimate < 0, 1,
  ifelse(BC_OE_bulk$FDR < 0.05 & BC_OE_bulk$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p3 <- EnhancedVolcano(BC_OE_bulk,
                      lab = BC_OE_bulk$Gene,
                      selectLab = c(	'Aox2', 'Cyp2f2', 'Unc119', 'Sord', 'Uqcrh', 'Srd5a2', 'Gchfr', 'Cyp2f2', 'Cyb5a', 'Ifit3', 'Isg15', 'Cd274', 'Olfr1049', 'Bst2', 'B2m', 'Ccl5', 'H2-Ab1', 'Tap1', 'Irf7', 'Mx2', 'Zbp1', 'Acod1', 'Junb', 'Cxcl10', 'Junb' ),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Olfactory Epithelium: Virus Absent vs. Virus Present',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-3, 3),
                      ylim = c(0, 2),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 20)
p3


#######################################################################################################

# OLFACTORY SENSORY NEURONS: MOCK VS. VIRUS ABSENT (Group A vs. Group B) 

keyvals.color <- ifelse(
  AB_OSN$Estimate < 0, '#AFAFB0',
  ifelse(AB_OSN$Estimate > 0, '#FB02FF', 'black'))

names(keyvals.color)[keyvals.color == '#FB02FF'] <- 'Enriched in Infected, virus absent'
names(keyvals.color)[keyvals.color == '#AFAFB0'] <- 'Enriched in Mock-infected'


keyvals.shape <- ifelse(
  AB_OSN$FDR <0.05 & AB_OSN$Estimate < 0, 1,
  ifelse(AB_OSN$FDR < 0.05 & AB_OSN$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p4 <- EnhancedVolcano(AB_OSN,
                      lab = AB_OSN$Gene,
                      selectLab = c('Olfr120', 'Bst2', 'Irf7', 'Lgals3bp', 'Aqr', 'B2m', 'Usp18', 'Ifitm3', 'Oasl2', 'Nasp', 'Tapbp', 'Cox8a', 'Ift80', 'Ndufa4', 'Rap1gds1', 'Nme4', 'Dalrd3', 'Ddx60', 'Irf7', 'Irf9'),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Olfactory Sensory Neurons: Mock vs. Virus Absent',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-2, 5),
                      ylim = c(0, 2),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 20)
p4


#######################################################################################################

# OLFACTORY SENSORY NEURONS: MOCK VS. VIRUS PRESENT (Group A vs. Group C)

keyvals.color <- ifelse(
  AC_OSN$Estimate < 0, '#AFAFB0',
  ifelse(AC_OSN$Estimate > 0, '#21FF06', 'black'))

names(keyvals.color)[keyvals.color == '#AFAFB0'] <- 'Enriched in Mock-infected'
names(keyvals.color)[keyvals.color == '#21FF06'] <- 'Enriched in Infected, virus present'


keyvals.shape <- ifelse(
  AC_OSN$FDR <0.05 & AC_OSN$Estimate < 0, 1,
  ifelse(AC_OSN$FDR < 0.05 & AC_OSN$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p5 <- EnhancedVolcano(AC_OSN,
                      lab = AC_OSN$Gene,
                      selectLab = c('Olfr59', 'Bst2', 'C1qb', 'B2m', 'Irf9', 'Rap1gds1', 'Isg15', 'Irf7', 'Oas1a', 'Adar', 'Tap1', 'Xaf1', 'Cd274', 'Ifi44', 'Tlr3', 'Bsg', 'Zbp1', 'H4c8', 'Cbx4', 'Hsp90aa1', 'Ank3', 'Rtp4', 'Ifitm3', "C1qb", "Cxcr4", 'Lgals3bp', 'Mpeg1'),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Olfactory Sensory Neurons: Mock vs. Virus Present',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-3, 8),
                      ylim = c(0, 4),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 10)
p5

#######################################################################################################

# OLFACTORY SENSORY NEURONS: VIRUS ABSENT VS. VIRUS PRESENT (Group B vs. Group C) 

keyvals.color <- ifelse(
  BC_OSN$Estimate < 0, '#FB02FF',
  ifelse(BC_OSN$Estimate > 0, '#21FF06', 'black'))

names(keyvals.color)[keyvals.color == '#FB02FF'] <- 'Enriched in Infected, virus absent'
names(keyvals.color)[keyvals.color == '#21FF06'] <- 'Enriched in Infected, virus present'


keyvals.shape <- ifelse(
  BC_OSN$FDR <0.05 & BC_OSN$Estimate < 0, 1,
  ifelse(BC_OSN$FDR < 0.05 & BC_OSN$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p6 <- EnhancedVolcano(BC_OSN,
                      lab = BC_OSN$Gene,
                      #selectLab = c(),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Olfactory Sensory Neurons: Virus Absent vs. Virus Present',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-3, 3),
                      ylim = c(0, 2),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 20)
p6

#######################################################################################################

# BULK GLOMERULI: MOCK VS. VIRUS ABSENT (Group A vs. Group B) 

keyvals.color <- ifelse(
  AB_GL_bulk$Estimate < 0, '#AFAFB0',
  ifelse(AB_GL_bulk$Estimate > 0, '#FB02FF', 'black'))

names(keyvals.color)[keyvals.color == '#AFAFB0'] <- 'Enriched in Mock-infected'
names(keyvals.color)[keyvals.color == '#FB02FF'] <- 'Enriched in Infected, virus absent'


keyvals.shape <- ifelse(
  AB_GL_bulk$FDR <0.05 & AB_GL_bulk$Estimate < 0, 1,
  ifelse(AB_GL_bulk$FDR < 0.05 & AB_GL_bulk$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p7 <- EnhancedVolcano(AB_GL_bulk,
                      lab = AB_GL_bulk$Gene,
                      selectLab = c('Olfr272', 'Hnrnpa2b1', 'Cck', 'Ncan', 'Ninj2', 'Hsp90b1', 'Irf7', 'Lcat', 'Gm1110', 'Irf9', 'Ifit3', 'Btbd19', 'Krt88', 'Bsg', 'Podxl'),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Glomeruli: Mock vs. Virus Absent',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-2, 2),
                      ylim = c(0, 2),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 20)
p7

#######################################################################################################

# SUSENTACULAR CELLS: MOCK vs. VIRUS ABSENT (Group A vs. Group B)

keyvals.color <- ifelse(
  AB_SUS$Estimate < 0, '#AFAFB0',
  ifelse(AB_SUS$Estimate > 0, '#FB02FF', 'black'))

names(keyvals.color)[keyvals.color == '#AFAFB0'] <- 'Enriched in Mock-infected'
names(keyvals.color)[keyvals.color == '#FB02FF'] <- 'Enriched in Infected, virus absent'


keyvals.shape <- ifelse(
  AB_SUS$FDR <0.05 & AB_SUS$Estimate < 0, 1,
  ifelse(AB_SUS$FDR < 0.05 & AB_SUS$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p8 <- EnhancedVolcano(AB_SUS,
                      lab = AB_SUS$Gene,
                      #selectLab = c(),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Susentacular cells: Mock vs. Virus Absent',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-2, 3),
                      ylim = c(0, 2),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 20)
p8

#######################################################################################################

# SUSENTACULAR CELLS: MOCK VS. VIRUS PRESENT (Group A vs. Group C) Susentacular cells

keyvals.color <- ifelse(
  AC_SUS$Estimate < 0, '#AFAFB0',
  ifelse(AC_SUS$Estimate > 0, '#21FF06', 'black'))

names(keyvals.color)[keyvals.color == '#AFAFB0'] <- 'Enriched in Mock-infected'
names(keyvals.color)[keyvals.color == '#21FF06'] <- 'Enriched in Infected, virus present'


keyvals.shape <- ifelse(
  AC_SUS$FDR <0.05 & AC_SUS$Estimate < 0, 1,
  ifelse(AC_SUS$FDR < 0.05 & AC_SUS$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p9 <- EnhancedVolcano(AC_SUS,
                      lab = AC_SUS$Gene,
                      #selectLab = c(),
                      x = "Estimate",
                      y = "FDR",
                      title = 'Susentacular cells: Mock vs. Virus Present',
                      titleLabSize = 14,
                      pCutoff = 0.05,
                      FCcutoff = 0.5,
                      pointSize = 4.0,
                      labSize = 4,
                      colCustom = keyvals.color,
                      shapeCustom = keyvals.shape,
                      shape = 1,
                      colAlpha = 0.9,
                      legendLabSize = (10), 
                      xlim = c(-2, 4),
                      ylim = c(0, 4),
                      axisLabSize = 15,
                      widthConnectors = 0.5,
                      drawConnectors = TRUE,
                      legendPosition = 'right',
                      subtitle = NULL,
                      caption = NULL,
                      gridlines.major = FALSE,
                      gridlines.minor = FALSE,
                      cutoffLineWidth = .1, 
                      legendIconSize = 3,
                      maxoverlapsConnectors = 10)
p9

#######################################################################################################

# SUSENTACULAR CELLS: VIRUS ABSENT VS. VIRUS PRESENT (Group B vs. Group C)

keyvals.color <- ifelse(
  BC_SUS$Estimate < 0, '#FB02FF',
  ifelse(BC_SUS$Estimate > 0, '#21FF06', 'black'))

names(keyvals.color)[keyvals.color == '#FB02FF'] <- 'Enriched in Infected, virus absent'
names(keyvals.color)[keyvals.color == '#21FF06'] <- 'Enriched in Infected, virus present'


keyvals.shape <- ifelse(
  BC_SUS$FDR <0.05 & BC_SUS$Estimate < 0, 1,
  ifelse(BC_SUS$FDR < 0.05 & BC_SUS$Estimate > 0, 2, 3))

names(keyvals.shape)[keyvals.shape == 1] <- 'Log2FC < 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 2] <- 'Log2FC > 0, FDR < 0.05'
names(keyvals.shape)[keyvals.shape == 3] <- 'NS'


p10 <- EnhancedVolcano(BC_SUS,
                       lab = BC_SUS$Gene,
                       #selectLab = c(),
                       x = "Estimate",
                       y = "FDR",
                       title = 'Susentacular cells: Virus Absent vs. Virus Present',
                       titleLabSize = 14,
                       pCutoff = 0.05,
                       FCcutoff = 0.5,
                       pointSize = 4.0,
                       labSize = 4,
                       colCustom = keyvals.color,
                       shapeCustom = keyvals.shape,
                       shape = 1,
                       colAlpha = 0.9,
                       legendLabSize = (10), 
                       xlim = c(-3, 3),
                       ylim = c(0, 2),
                       axisLabSize = 15,
                       widthConnectors = 0.5,
                       drawConnectors = TRUE,
                       legendPosition = 'right',
                       subtitle = NULL,
                       caption = NULL,
                       gridlines.major = FALSE,
                       gridlines.minor = FALSE,
                       cutoffLineWidth = .1, 
                       legendIconSize = 3,
                       maxoverlapsConnectors = 20)
p10

#######################################################################################################

# EXPORT HI-RES PLOTS

## Olfactory epithelium
#ggsave(file = "OE_AB.tiff", plot = p1, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )
#ggsave(file = "OE_AC.tiff", plot = p2, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )
#ggsave(file = "OE_BC.tiff", plot = p3, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )

# Olfactory sensory neurons
#ggsave(file = "OSN_AB.tiff", plot = p4, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )
#ggsave(file = "OSN_AC.tiff", plot = p5, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )
#ggsave(file = "OSN_BC.tiff", plot = p6, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )

# Glomeruli
#ggsave(file = "GL_AB.tiff", plot = p7, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )

# Susentacular cells
#ggsave(file = "SUS_AB.tiff", plot = p8, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )
#ggsave(file = "SUS_AC.tiff", plot = p9, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )
#ggsave(file = "SUS_BC.tiff", plot = p10, device = "tiff", path = "./write/figures/", width = 8, height = 4, dpi = 600 )

#######################################################################################################