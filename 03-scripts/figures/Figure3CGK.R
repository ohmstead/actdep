## ----Fig3CGK
# Script to plot heatmaps from figure 1, with gene columns organized by tau.
library(ComplexHeatmap)

source('03-scripts/R/seq_functions.R')

gigaclasses <- LoadSubclassesToUse(nuclei, as_gigaclasses = T)
activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# establish subclasses to use ----------------------------------------

ieg_symbols <- c(LoadGeneList('IEG'), 'Rn7sk', 'Midn', 
                 'Etv1', 'Rcan2', 'Hs3st2',
                 'Tac1', 'Socs2', 'Daam2', 'Akap5', 'Frmd6', 'Crhbp', 'Gm49673',
                 'Map3k19', 'Usp53', 'Slco1c1', 'Kcnn2')

# loop through cell superclasses
plots_tau <- list()
for (gigaclass_str in names(gigaclasses)) {
  subclass_list <- gigaclasses[[gigaclass_str]]
  
  # load df_expression
  load_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
  fpath <- glue("{load_dir}/0_df_expression_{gigaclass_str}.csv")
  df_expression <- read_csv(fpath, show_col_types = F)
  
  # calculate tau
  df_sub <- df_expression |>
    filter(activity_condition == 'EE30m',
           classification == 'ERG',
           direction == 'up'
           ) |>
    group_by(gene) |> 
    mutate(tau = tau(log2FC_calculated, direction = 'up')) |>
    select(gene, subclass, log2FC_calculated, tau)
  
  gene_tau <- df_sub |> 
    group_by(gene) |> 
    slice_head() |> 
    select(gene, tau) |> 
    deframe()
  
  # Pivot to wide format so that rows are subclasses and columns are genes.
  mat <- df_sub |>
    select(-tau) |> 
    pivot_wider(names_from = gene, values_from = log2FC_calculated) |>
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
  mat_ordered <- mat[, gene_order, drop = FALSE]  # set gene cols
  
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
               col = circlize::colorRamp2(c(-3, 0, 3), hcl_palette = 'RdBu', reverse = T),
               row_order = subclass_list,
               top_annotation = top_anno,
               bottom_annotation = bottom_anno,
               height = unit(7.5*length(subclass_list), 'mm'),
               cluster_rows = F,
               cluster_columns = F,
               show_row_names = T,
               show_column_names = F
  )
  draw(p, heatmap_legend_side = "bottom")
  
  plots_tau[[paste(gigaclass_str, 'ERG', sep = "_")]] <- p
}

# InteractiveComplexHeatmap::htShiny(plots_tau$glia_ERG, output_ui_float = T)


# save ----------------------------------------
if (exists('SAVE_PLOTS') & SAVE_PLOTS) {
  save_dir <- "05-results/Figure3/raw_R_plots"
  for (current_plot in names(plots_tau)) {
    p <- plots_tau[[current_plot]]
    plot_path <- glue("{save_dir}/PNG_{current_plot}.png")
    
    # png
    png(plot_path, width = 8, height = 5, units = "in", res = 900)
    draw(p, heatmap_legend_side = 'bottom')
    dev.off()
    
    # svg
    svgsave(plot = p,
            filename = glue("SVG_{current_plot}.svg"), 
            path = save_dir, 
            width = 8, height = 5)
  }
}
## ----