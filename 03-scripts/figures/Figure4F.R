## ----Fig4F
# This script produces correlation matricies for IEGs and other DEGs in CA1 supertypes
source('03-scripts/R/seq_functions.R')

library(purrr)
library(circlize)
library(ComplexHeatmap)

if(!exists('nuclei')) {nuclei <- LoadDataset("Dec2024")}

# fxns --------------------------------------------
MakeHMobject <- function(mat, plot_title) {
  hm <- Heatmap(
    mat,
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
    rect_gp = gpar(type = "none"), 
    column_dend_side = "top",
    show_row_dend = F,
    show_heatmap_legend = F,
    height = unit(5, "in"),
    width  = unit(5, "in"),
    column_dend_height = unit(0.2, "in"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      if(as.numeric(x) > 1 - as.numeric(y) - 1e-6) {
        grid.rect(x, y, w, h, gp = gpar(fill = fill, col = 'black'))
      }
    }
  )
  # hm <- Heatmap(
  #   mat,
  #   col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  #   column_title = plot_title,
  #   cluster_columns = T,
  #   cluster_rows = T,
  #   show_heatmap_legend = F,
  #   width = 900,
  #   height = 900,
  #   row_names_gp = gpar(fontsize = 20),
  #   column_names_gp = gpar(fontsize = 20),
  # )
  return(hm)
}

save_path <- glue("05-results/Figure4/raw_R_plots/supertype_coexpression")
SaveHMobject <- function(hm, save_path, plot_title) {
  if (!exists('SAVE_PLOTS')) {return()} else { if (SAVE_PLOTS == F) {return()}}
  # png
  png(glue("{save_path}/{plot_title}.png"), 
      width = 9, height = 10, units = 'in', res = 900, bg = 'white')
  draw(hm)
  dev.off()
  # svg
  svgsave(plot = hm,
          savedir = save_path,
          filename = plot_title,
          w = 9, h = 10)
}


# get supertypes --------------------------------------------
IEG_symbols <- LoadGeneList("IEG")

if (!exists('FindActiveCells_CA1_outputs')) {
  FindActiveCells_CA1_outputs <- FindActiveCells(
    nuclei, 
    subclass = '016 CA1-ProS Glut', 
    gene_list = IEG_symbols, 
    gene_threshold = 3
  )
}
df_active_cells_IEGs <- FindActiveCells_CA1_outputs$df_active_cells |> 
  select(-activity_condition)  # redundant col for the merge
  

# plot supertypes --------------------------------------------
supertypes <- c(
  "0069 CA1-ProS Glut_1",
  "0070 CA1-ProS Glut_2",
  "0072 CA1-ProS Glut_4"
)

# make list of seurat objects for each supertype
supertype_seurats <- supertypes |> 
  set_names() |> 
  map(~ subset(nuclei, supertype_name == .x), .progress = F)

for (supertype in supertypes) {
  nuclei_supertype <- supertype_seurats[[supertype]]
  
  # let's get the active cells
  meta <- nuclei_supertype@meta.data |> 
    left_join(df_active_cells_IEGs, by = c("barcode" = "cell")) |> 
    relocate(barcode)
  
  nuclei_supertype <- AddMetaData(nuclei_supertype, meta)
  
  active_barcodes <- df_active_cells_IEGs |> 
    filter(active_binary == TRUE) |>
    distinct(cell) |> 
    pull(cell)
  
  active_cells <- nuclei_supertype |> subset(barcode %in% active_barcodes)
  
  # ACTIVE cells
  subset_cells_EE30m <- active_cells |> subset(activity_condition == 'EE30m')
  mat_coexpression <- GetCorrData(subset_cells_EE30m,
                                  specific_condition = 'EE30m',
                                  gene_list = IEG_symbols)
  
  plot_title <- glue("{supertype}__EE30m_active_cells") |> 
    str_replace_all(' ', '_')
  hm <- MakeHMobject(mat_coexpression, plot_title)
  draw(hm)

  SaveHMobject(hm, save_path, plot_title)
}
## ----