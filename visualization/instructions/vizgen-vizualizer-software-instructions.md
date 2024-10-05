# MERSCOPE “Vizualizer” Software Quick Guide

## Installing MERSCOPE Vizualizer
1. Log-in to Vizgen’s online portal.
2. Download and install “MERSCOPE Vizualizer” for Windows / macOS: [MERSCOPE Vizualizer](https://portal.vizgen.com/resources/software).

## Importing Data
1. Open MERSCOPE Vizualizer software.
2. **File > Assemble Experiment**
   - Enter experiment name, select output folder.
   - Select Datasets: import `.vzg` file associated with experiment.
   - Click Build.
3. **File > Import Data > AnnData**
   - Select custom `.hdf5` file.
   - Select annotations under “CATEGORY” (i.e., “leiden”, “cell_type”).
   - Click Import.

## Viewing Cells / Clusters
- Under “Cells/Spatial Features” (top left panel of GUI), change the “Active category” to your desired annotation mode:
  - **Leiden**: Cell clusters based on the Leiden clustering algorithm (successor to “Louvain” clustering).
  - **Cell_type**: Manual Leiden cluster assignment based on gene expression.
  - Additional cell groupings/annotations may appear in this section.
- When viewing cells, uncheck “Unclustered” to remove cells/observations that did not pass QC metrics and filtering.
- Ensure “Cells” box is checked in the center panel above the interactive image of the tissue slice.

## Viewing Genes
- Transcript information is located under “Genes,” below the “Spatial Features” section.
- Spatial information for each gene is available in this section.
- Ensure “Transcripts” box is checked in the center panel above the interactive image of the tissue slice.

## Drawing custom ROI
- Use polygon tool to draw custom ROI and export geometries to CSV for compatibility with analysis in Python.