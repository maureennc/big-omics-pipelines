library(readr)
library(EnhancedVolcano)
library(ggplot2)

##############################################################################

# IMPORT DATA

## Define paths
root_dir <- '/Users/maureen/Documents/projects/lukens-lab/ana/2024_tbi-snrna-seq'
data_dir <- file.path(root_dir, 'analysis', '1', 'results', 'spreadsheets', 'de', 'cell_type', 'filtered')

## Read and prepare data
read_and_prepare_data <- function(comparison, fc_cutoff) {
  file_path <- file.path(data_dir, paste(comparison, "comparison-cell_type-filtered.csv", sep="-"))
  data <- read_csv(file_path)
  cell_types <- unique(data$cell_type)
  list_data <- list()
  
  for (cell_type in cell_types) {
    sub_data <- data[data$cell_type == cell_type, ]
    filtered_data <- sub_data[sub_data$qval < 0.05 & abs(sub_data$log2fc) >= fc_cutoff, ]
    list_data[[cell_type]] <- list(data = sub_data, filtered_genes = filtered_data$gene)
  }
  
  return(list_data)
}

##############################################################################

# FUNCTION TO GENERATE VOLCANO PLOTS

generate_volcano_plots <- function(data_list, comparison, fc_cutoff) {
  plots <- list()
  
  for (cell_type in names(data_list)) {
    data <- data_list[[cell_type]]$data
    filtered_genes <- data_list[[cell_type]]$filtered_genes
    plot <- EnhancedVolcano(data,
                            lab = data$gene,
                            x = 'log2fc',
                            y = 'pval',
                            pCutoffCol = 'qval',
                            title = paste(cell_type, comparison, sep=" - "),
                            pCutoff = 0.05,
                            FCcutoff = fc_cutoff,
                            pointSize = 3.0,
                            labSize = 7,
                            selectLab = filtered_genes,
                            drawConnectors = TRUE,
                            widthConnectors = 0.75,
                            legendLabSize = 10,
                            max.overlaps = 20)
    plots[[cell_type]] <- plot
  }
  
  return(plots)
}

##############################################################################

# DEFINE PARAMETERS

comparisons <- c('AC', 'CE', 'AB', 'AD', 'CD', 'DF')
fc_cutoffs <- c(0.5, 1)

##############################################################################

# CREATE AND SAVE PLOTS

for (comparison in comparisons) {
  for (fc_cutoff in fc_cutoffs) {
    data_list <- read_and_prepare_data(comparison, fc_cutoff)
    plots <- generate_volcano_plots(data_list, comparison, fc_cutoff)
    
    ## Save plots
    save_dir <- file.path(root_dir, 'analysis', '1', 'results', 'figures', 'volcano', paste("logfc", fc_cutoff, comparison, sep="_"))
    if (!dir.exists(save_dir)) {
      dir.create(save_dir, recursive = TRUE)
    }
    
    for (cell_type in names(plots)) {
      plot_file <- file.path(save_dir, paste(cell_type, comparison, "fc", fc_cutoff, "volcano.pdf", sep="_"))
      ggsave(plot_file, plots[[cell_type]], width = 6, height = 6, dpi = 300)
    }
  }
}

##############################################################################