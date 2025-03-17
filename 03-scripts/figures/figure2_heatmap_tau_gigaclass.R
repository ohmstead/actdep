library(ggplot2)
library(dplyr)
library(tidyr)

library(Seurat)
library(ComplexHeatmap)

source('03-scripts/R/seq_functions.R')
activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# establish subclasses to use ----
print('Subsetting Seurat object...')
subclass_sets <- LoadSubclassesToUse(nuclei, as_gigaclasses = T)
seurat_subsets <- LoadDataset("Dec2024", as_gigaclasses = T)

ieg_symbols <- c(LoadGeneList('IEG'), 'Rn7sk', 'Rcan2', 'Tac1', 'Map3k19')

# loop through cell superclasses
plots_tau <- list()
for (set_name in names(seurat_subsets)) {
  subclass_list <- subclass_sets[[set_name]]
  nuclei_subclass <- seurat_subsets[[set_name]]
  nuclei_subclass$subclass_name |> table()
  df_expression <- read_csv(glue("04-analysis/df_expression/df_expression.seurat/df_expression_{set_name}.csv")) |> 
    print()
  
  df_sub <- df_expression |>
    filter(activity_condition == 'EE30m',
           classification == 'ERG',
           direction == 'up'
           ) |>
    group_by(gene) |> 
    mutate(tau = tau(log2FoldChange.shrink, direction = 'up')) |>
    select(gene, subclass, log2FoldChange.shrink, tau)
  
  gene_tau <- df_sub |> 
    group_by(gene) |> 
    slice_head() |> 
    select(gene, tau) |> 
    deframe()
  
  # Pivot to wide format so that rows are subclasses and columns are genes.
  mat <- df_sub |>
    select(-tau) |> 
    pivot_wider(names_from = gene, values_from = log2FoldChange.shrink) |>
    column_to_rownames("subclass") |>
    as.matrix()
  
  # # Compute tau for each gene across the subclasses (columns of the matrix)
  # gene_tau <- apply(mat, 2, function(x) tau(x, direction = "up"))
  # 
  # # Force tau between 0 and 1 using indexing
  # gene_tau[gene_tau < 0] <- 0
  # gene_tau[gene_tau > 1] <- 1
  
  # Order genes by tau (highest to lowest)
  gene_order <- names(sort(gene_tau, decreasing = TRUE))
  bar_colors <- ifelse(gene_tau[gene_order] > 0.8, "red", "gray60")
  mat_ordered <- mat[, gene_order, drop = FALSE]
  
  # (Optional) Create a top annotation showing tau as a barplot
  top_anno <- HeatmapAnnotation(
    tau = anno_barplot(gene_tau[gene_order], 
                       border = FALSE, 
                       bar_width = 0.8,
                       height = unit(1.5, 'cm'),
                       gp = gpar(fill = bar_colors))
  )
  
  bottom_anno <- HeatmapAnnotation(
    iegs = anno_mark(at = which(gene_order %in% ieg_symbols), 
                     labels = intersect(colnames(mat_ordered), ieg_symbols), 
                     side = "bottom"),
    show_annotation_name = FALSE
  )
  
  # plot ----------------------------------------------------------------
  p <- Heatmap(mat_ordered,
               name = "log2FC",
               col = circlize::colorRamp2(c(-2, 0, 2), hcl_palette = 'Blue-Red 2'),
               top_annotation = top_anno,
               bottom_annotation = bottom_anno,
               cluster_rows = FALSE,
               cluster_columns = FALSE,
               show_row_names = FALSE,
               show_column_names = FALSE,
               # column_title = paste(set_name, 'ERG', "tau-organized heatmap"),
               heatmap_legend_param = list(title = "log2FC"))
  draw(p, heatmap_legend_side = "bottom")
  
  plots_tau[[paste(set_name, 'ERG', sep = "_")]] <- p
}

htShiny(plots_tau[[1]])