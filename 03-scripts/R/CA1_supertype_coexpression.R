# for each supertype, I want to create their own complex heatmap
library(tidyverse)
library(Seurat)
library(glue)

# load in nuclei
# nuclei <- readRDS("04-analysis/Seurats/May2024/seurat.Rds")
nuclei <- nuclei |> subset(sublibrary == 1)

# first, get all unique supertypes
supertypes <- nuclei@meta.data |> 
  filter(subclass_name == '016 CA1-ProS Glut') |> 
  distinct(supertype_name) |> 
  pull(supertype_name)

conditions <- nuclei@meta.data |> 
  distinct(activity_condition) |> 
  pull(activity_condition)

IEG_symbols <- LoadGeneList("IEG")


# supertypes ACTIVE/INACTIVE cells ----
for (supertype in supertypes) {
  cells <- nuclei |> 
    subset(supertype_name == supertype)
  
  # let's find the active cells
  output_vars <- FindActiveCells(cells, gene_list = IEG_symbols, gene_threshold = 3)
  df_active_cells <- output_vars$df_active_cells |> 
    select(-activity_condition)  # redundant col for the merge
  
  meta <- cells@meta.data |> 
    left_join(output_vars$df_active_cells, by = c("bc_wells" = "cell")) |> 
    relocate(bc_wells)
  
  cells <- AddMetaData(cells, meta)
  
  active_barcodes <- output_vars$df_active_cells |> 
    filter(active_binary == TRUE) |>
    distinct(cell) |> 
    pull(cell)
  
  active_cells <- cells |> 
    subset(bc_wells %in% active_barcodes)
  inactive_cells <- cells |> 
    subset(bc_wells %in% active_barcodes, invert = TRUE)
  
  for (condition in conditions) {
    tryCatch(
      {
        # do it for active cells
        subset_cells <- active_cells |> 
          subset(activity_condition == condition)
        
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{supertype}__{condition}_active_cells")
        save_path <- glue("04-analysis/coexpression_heatmaps_supertype/")
        PlotComplexHeatmap(mat_coexpression, plot_title = plot_title, save_path = save_path)
        
        # do it for inactive cells
        subset_cells <- inactive_cells |> 
          subset(activity_condition == condition)
        
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{supertype}__{condition}_inactive_cells")
        save_path <- glue("04-analysis/coexpression_heatmaps_supertype/")
        PlotComplexHeatmap(mat_coexpression, plot_title = plot_title, save_path = save_path)
      },
      error = function(e) {
        print(glue("{supertype}, condition {condition} has no active cells"))
      }
    )
  }
}


# supertypes ALL cells ----
for (supertype in supertypes) {
  cells <- nuclei |> 
    subset(supertype_name == supertype)
  
  for (condition in conditions) {
    tryCatch(
      {
        subset_cells <- cells |> 
          subset(activity_condition == condition)
        
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{supertype}__{condition}_all_cells")
        save_path <- glue("04-analysis/coexpression_heatmaps_supertype/")
        PlotComplexHeatmap(mat_coexpression, plot_title = plot_title, save_path = save_path)
      },
      error = function(e) {
        print(glue("{supertype}, condition {condition} doesn't have enough cells"))
      }
    )
  }
}


# clusters ACTIVE/INACTIVE cells ----
clusters <- nuclei@meta.data |> 
  filter(subclass_name == '016 CA1-ProS Glut') |> 
  filter(str_detect(cluster_name, "262|263")) |> 
  distinct(cluster_name) |> 
  arrange(cluster_name) |>
  pull(cluster_name)

for (cluster in clusters) {
  cells <- nuclei |> 
    subset(cluster_name == cluster)
  
  # let's find the active cells
  output_vars <- FindActiveCells(cells, gene_list = IEG_symbols, gene_threshold = 3)
  df_active_cells <- output_vars$df_active_cells |> 
    select(-activity_condition)  # redundant col for the merge
  
  meta <- cells@meta.data |> 
    left_join(output_vars$df_active_cells, by = c("bc_wells" = "cell")) |> 
    relocate(bc_wells)
  
  cells <- AddMetaData(cells, meta)
  
  active_barcodes <- output_vars$df_active_cells |> 
    filter(active_binary == TRUE) |>
    distinct(cell) |> 
    pull(cell)
  
  active_cells <- cells |> 
    subset(bc_wells %in% active_barcodes)
  inactive_cells <- cells |> 
    subset(bc_wells %in% active_barcodes, invert = TRUE)
  
  for (condition in conditions) {
    tryCatch(
      {
        # do it for active cells
        subset_cells <- active_cells |> 
          subset(activity_condition == condition)
        
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{cluster}__{condition}_active_cells")
        save_path <- glue("04-analysis/coexpression_heatmaps_cluster/")
        PlotComplexHeatmap(mat_coexpression, plot_title = plot_title, save_path = save_path)
        
        # do it for inactive cells
        subset_cells <- inactive_cells |> 
          subset(activity_condition == condition)
        
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{cluster}__{condition}_inactive_cells")
        save_path <- glue("04-analysis/coexpression_heatmaps_cluster/")
        PlotComplexHeatmap(mat_coexpression, plot_title = plot_title, save_path = save_path)
      },
      error = function(e) {
        print(glue("{cluster}, condition {condition} has no inactive cells"))
      }
    )
  }
}


# clusters ALL cells ----
for (cluster in clusters) {
  cells <- nuclei |> 
    subset(cluster_name == cluster)
  
  for (condition in conditions) {
    tryCatch(
      {
        subset_cells <- cells |> 
          subset(activity_condition == condition)
        
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{cluster}__{condition}_all_cells")
        save_path <- glue("04-analysis/coexpression_heatmaps_cluster/")
        PlotComplexHeatmap(mat_coexpression, plot_title = plot_title, save_path = save_path)
      },
      error = function(e) {
        print(glue("{cluster}, condition {condition} doesn't have enough cells"))
      }
    )
  }
}



tmp <- nuclei |>
  subset(supertype_name == '0069 CA1-ProS Glut_1')

res <- FindActiveCells(tmp, gene_list = IEG_symbols, gene_threshold = 3)

df_active_cells <- res$df_active_cells |>
  select(-activity_condition) |>   # redundant col for the merge
  filter(num_upregd_genes == 4)  # only keep cells with at least 1 gene upregulated)

active_barcodes <- df_active_cells |>
  filter(active_binary == TRUE) |>
  distinct(cell) |>
  pull(cell)
tmp.subset <- tmp |>
  subset(bc_wells %in% active_barcodes)
tmp.mat <- GetCorrData(tmp.subset, specific_condition = 'SE', gene_list = IEG_symbols)
PlotComplexHeatmap(tmp.mat, plot_title = '0069 CA1-ProS Glut_1__SE_4ieg_active_cells', save_path = '~/Downloads/')
