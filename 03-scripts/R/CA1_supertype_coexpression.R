# This script produces correlation matricies for IEGs and other DEGs in CA1 supertypes

library(tidyverse)
library(glue)

library(Seurat)
library(ComplexHeatmap)

# load in nuclei
nuclei <- LoadDataset("Dec2024")

# first, get all unique supertypes
supertypes <- nuclei@meta.data |> 
  filter(subclass_name == '016 CA1-ProS Glut') |> 
  distinct(supertype_name) |> 
  pull(supertype_name)

conditions <- nuclei@meta.data |> 
  distinct(activity_condition) |> 
  pull(activity_condition)

IEG_symbols <- LoadGeneList("IEG")


# fxns: HM plotting and saving ------------------------------------------------
MakeHMobject <- function(mat, plot_title) {
  hm <- Heatmap(
    mat,
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    column_title = plot_title,
    cluster_columns = T,
    cluster_rows = T,
    show_heatmap_legend = F,
    width = 900,
    height = 900,
    row_names_gp = gpar(fontsize = 20),
    column_names_gp = gpar(fontsize = 20),
  )
  return(hm)
}

SaveHMobject <- function(hm, save_path, plot_title) {
  file_path <- glue("{save_path}/{plot_title}.svg")
  
  svg(file_path, width = 10, height = 9)
  draw(hm)
  dev.off()
}


# supertypes ACTIVE/INACTIVE cells ------------------------------------------------
for (supertype in supertypes) {
  nuclei_supertype <- nuclei |> subset(supertype_name == supertype)
  save_path <- glue("04-analysis/coexpression_heatmaps_supertype/")
  
  # let's find the active cells
  output_vars <- FindActiveCells(nuclei_supertype, gene_list = IEG_symbols, gene_threshold = 3)
  df_active_cells <- output_vars$df_active_cells |> 
    select(-activity_condition)  # redundant col for the merge
  
  meta <- nuclei_supertype@meta.data |> 
    left_join(output_vars$df_active_cells, by = c("barcode" = "cell")) |> 
    relocate(barcode)
  
  nuclei_supertype <- AddMetaData(nuclei_supertype, meta)
  
  active_barcodes <- output_vars$df_active_cells |> 
    filter(active_binary == TRUE) |>
    distinct(cell) |> 
    pull(cell)
  
  active_cells <- nuclei_supertype |> subset(barcode %in% active_barcodes)
  inactive_cells <- nuclei_supertype |> subset(barcode %in% active_barcodes, invert = TRUE)
  
  for (condition in conditions) {
    tryCatch(
      {
        # ACTIVE cells
        subset_cells <- active_cells |> subset(activity_condition == condition)
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{supertype}__{condition}_active_cells") |> str_replace_all(' ', '_')
        hm <- MakeHMobject(mat_coexpression, plot_title)
        draw(hm)
        SaveHMobject(hm, save_path, plot_title)
        
        # INACTIVE cells
        subset_cells <- inactive_cells |> subset(activity_condition == condition)
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{supertype}__{condition}_inactive_cells") |> str_replace_all(' ', '_')
        hm <- MakeHMobject(mat_coexpression, plot_title)
        draw(hm)
        SaveHMobject(hm, save_path, plot_title)
      },
      error = function(e) {
        print(glue("{supertype}, condition {condition} has no active cells"))
      }
    )
  }
}


# supertypes ALL cells ------------------------------------------------
for (supertype in supertypes) {
  nuclei_supertype <- nuclei |> subset(supertype_name == supertype)
  save_path <- glue("04-analysis/coexpression_heatmaps_supertype")
  
  for (condition in conditions) {
    tryCatch(
      {
        subset_cells <- nuclei_supertype |> subset(activity_condition == condition)
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{supertype}__{condition}_all_cells") |> str_replace_all(' ', '_')
        hm <- MakeHMobject(mat_coexpression, plot_title)
        draw(hm)
        SaveHMobject(hm, save_path, plot_title)
      },
      error = function(e) {
        print(glue("{supertype}, condition {condition} doesn't have enough cells"))
      }
    )
  }
}


# clusters ACTIVE/INACTIVE cells ------------------------------------------------
clusters <- nuclei@meta.data |> 
  filter(subclass_name == '016 CA1-ProS Glut') |> 
  filter(str_detect(cluster_name, "262|263")) |> 
  distinct(cluster_name) |> 
  arrange(cluster_name) |>
  pull(cluster_name)

for (cluster in clusters) {
  nuclei_supertype <- nuclei |> subset(cluster_name == cluster)
  save_path <- glue("04-analysis/coexpression_heatmaps_cluster")
  
  # let's find the active cells
  output_vars <- FindActiveCells(nuclei_supertype, gene_list = IEG_symbols, gene_threshold = 3)
  df_active_cells <- output_vars$df_active_cells |> 
    select(-activity_condition)  # redundant col for the merge
  
  meta <- nuclei_supertype@meta.data |> 
    left_join(output_vars$df_active_cells, by = c("bc_wells" = "cell")) |> 
    relocate(bc_wells)
  
  nuclei_supertype <- AddMetaData(nuclei_supertype, meta)
  
  active_barcodes <- output_vars$df_active_cells |> 
    filter(active_binary == TRUE) |>
    distinct(cell) |> 
    pull(cell)
  
  active_cells <- nuclei_supertype |> subset(bc_wells %in% active_barcodes)
  inactive_cells <- nuclei_supertype |> subset(bc_wells %in% active_barcodes, invert = TRUE)
  
  for (condition in conditions) {
    tryCatch(
      {
        # ACTIVE cells
        subset_cells <- active_cells |> subset(activity_condition == condition)
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{cluster}__{condition}_active_cells") |> str_replace_all(' ', '_')
        hm <- MakeHMobject(mat_coexpression, plot_title)
        draw(hm)
        SaveHMobject(hm, save_path, plot_title)
        
        # INACTIVE cells
        subset_cells <- inactive_cells |> subset(activity_condition == condition)
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{cluster}__{condition}_inactive_cells") |> str_replace_all(' ', '_')
        hm <- MakeHMobject(mat_coexpression, plot_title)
        draw(hm)
        SaveHMobject(hm, save_path, plot_title)
      },
      error = function(e) {
        print(glue("{cluster}, condition {condition} has no inactive cells"))
      }
    )
  }
}


# clusters ALL cells ------------------------------------------------
for (cluster in clusters) {
  nuclei_supertype <- nuclei |> subset(cluster_name == cluster)
  save_path <- glue("04-analysis/coexpression_heatmaps_cluster")
  
  for (condition in conditions) {
    tryCatch(
      {
        subset_cells <- nuclei_supertype |> 
          subset(activity_condition == condition)
        
        mat_coexpression <- GetCorrData(subset_cells,
                                        specific_condition = condition,
                                        gene_list = IEG_symbols)
        
        plot_title <- glue("{cluster}__{condition}_all_cells") |> str_replace_all(' ', '_')
        
        hm <- MakeHMobject(mat_coexpression, plot_title)
        draw(hm)
        SaveHMobject(hm, save_path, plot_title)
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

hm <- MakeHMobject(mat_coexpression, plot_title)
draw(hm)