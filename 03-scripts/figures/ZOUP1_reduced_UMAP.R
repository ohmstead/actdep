library(ggplot2)
library(dplyr)
library(readr)
library(glue)
library(patchwork)

library(Seurat)

source("03-scripts/R/seq_functions.R")
activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# establish subclasses to use ----------------------------------------------------
print('Subsetting Seurat object...')
gigaclasses <- LoadSubclassesToUse(nuclei, as_gigaclasses = T)
seurat_subsets <- LoadDataset("Dec2024", as_gigaclasses = T)


# load DEGs ---------------------------------------------------------------------
print("Getting DEGs for all subclasses and contrasts...")
gene_list <- LoadGeneList('DEGs')
gene_list <- LoadGeneList("IEG")

# re-embed cells in this smaller gene space -----------------------------------
for (gigaclass in names(gigaclasses)) {
  print(glue("Creating subsetted Seurat object for {gigaclass}..."))
  print(glue("Re-embedding {gigaclass} cells in DEG space..."))
  nuclei_tmp_processed <- seurat_subsets[[gigaclass]] |> 
    subset(activity_condition %in% c('SE', 'EE30m', 'EE6h')) |> 
    RunPCA(features = gene_list) |> 
    FindNeighbors(dims = 1:9) |> 
    FindClusters(resolution = 0.8, random.seed = 17) |> 
    RunUMAP(dims = 1:9, seed.use = 17)
  
  
  # plot reduced UMAP ------------------------------------------------------------
  p1 <- DimPlot(nuclei_tmp_processed, reduction = 'umap', group.by = 'subclass_name') +
    scale_color_manual(values = subclass_colors) +
    labs(title = glue("{gigaclass} re-embed")) +
    theme_void() +
    theme(legend.position = 'none')
  
  p2 <- DimPlot(nuclei_tmp_processed, reduction = 'umap', group.by = 'activity_condition') +
    scale_color_manual(values = activity_colors) +
    labs(title = glue("{gigaclass} re-embed")) +
    theme_void() +
    theme(legend.position = 'none')
  
  p3 <- DimPlot(nuclei_tmp_processed, reduction = 'umap') +
    labs(title = glue("{gigaclass} re-embed")) +
    theme_void()
  
  print(p1 + p2 + p3)
  
  plots <- list()
  for (gene in gene_list) {
    print(glue("Plotting {gene}..."))
    
    p1 <- FeaturePlot(nuclei_tmp_processed, features = gene) +
      labs(title = glue("{gigaclass} {gene}")) +
      theme_void() +
      theme(legend.position = 'none')
    
    plots[[length(plots) + 1]] <- p1
  }
  p_mega <- wrap_plots(plots, ncol = 3, nrow = 5)
  print(p_mega)
}
