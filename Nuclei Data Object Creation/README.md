# Generating the sNucRNA-Seq Seurat Object

1. Download the CellRanger output files from GEO for each of the four regions: AD, AV, PD, PV
2. Run each region's Rmd notebook to generate an initial single-region Seurat object
3. Run the Integration notebook to combine the four regions together into the main Seurat object, with which we find clusters & generate the main UMAP plot.
