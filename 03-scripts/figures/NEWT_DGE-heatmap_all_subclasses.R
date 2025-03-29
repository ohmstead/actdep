# load libs and data ----
source('03-scripts/R/seq_functions.R')

library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(readxl)

activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# establish subclasses to use ----------------------------------------------------
print('Subsetting Seurat object...')
gigaclasses <- LoadSubclassesToUse(nuclei, as_gigaclasses = T)
# seurat_gigaclasses <- LoadDataset('Dec2024', as_gigaclass = T)

# loop thru sets of subclasses
plots <- list()
list_df_expression <- list(excitatory = NULL, inhibitory = NULL, glia = NULL)
for (gigaclass_str in names(gigaclasses)) {
  subclass_list <- gigaclasses[[gigaclass_str]]
  
  
  # read in df_expression ----------------------------------------------------
  print("Loading ERGs/LRGs (defined in classify_ERG-LRG.R)...")
  load_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
  
  df_expression <- read_csv(glue("{load_dir}/0_df_expression_{gigaclass_str}.csv"),show_col_types = F)
  df_classification <- read_csv(glue("{load_dir}/0_DEG_classifications_{gigaclass_str}.csv"),show_col_types = F)
  df_distinct_DEGs <- df_classification
  
  # get TF genes ----------------------------------------------------------------------------
  df_gene_levels <- df_classification
  
  # determine whether genes are TFs
  mouse_TFs <- read_excel("02-data/published_data/Mus_musculus_TF.xlsx") |> 
    select(Symbol) |> 
    pull()
  
  df_gene_levels <- df_gene_levels |>
    mutate(TF = ifelse(gene %in% mouse_TFs, "1", "0"))
  
  # combine classification and TF for anno purposes
  df_gene_levels$classification_TF <- paste(df_gene_levels$classification, df_gene_levels$TF, sep = "_")
  
  # re-level variables for plotting ----
  print("Re-leveling variables for plotting...")
  expression_matrix <- df_expression |> 
    filter(activity_condition != 'SE') |> 
    mutate(subclass = fct(subclass, levels = subclass_list)) |> 
    arrange(gene, classification, subclass) |> 
    select(gene, subclass_by_activity_condition, log2FC_calculated) |> 
    pivot_wider(names_from = 'gene', values_from = 'log2FC_calculated') |> 
    column_to_rownames('subclass_by_activity_condition') |> 
    as.data.frame() |> 
    as.matrix()
  
  expression_matrix <- expression_matrix[,df_gene_levels$gene]
  mat = list(
    ERG = expression_matrix[
      str_detect(rownames(expression_matrix), 'EE30m'),
      df_gene_levels$classification == 'ERG'
      ],
    LRG = expression_matrix[
      str_detect(rownames(expression_matrix), 'EE6h'),
      df_gene_levels$classification == 'LRG'
      ]
  )
  
  
  # make annotation objects ----------------------------------------------------------------------------
  print('Making annotation objects and plotting')
  df_looping <- data.frame(
    gene_classification = c('ERG', 'LRG'),
    condition           = c('EE30m', 'EE6h')
    # seriation           = I(list(erg_seriation, lrg_seriation))
    )
  for (k in 1:2) { # ERG/LRG
    looping <- df_looping[k,]
    plotting_matrix <- mat[[k]]
    
    df_gene_annotation <- df_gene_levels |> 
      filter(classification == looping$gene_classification) |> 
      mutate(gene = factor(gene, levels = colnames(plotting_matrix))) |> 
      arrange(gene)
      
    df_subclass_annotation <- df_expression |>
      filter(classification == looping$gene_classification) |>
      filter(activity_condition == looping$condition) |> 
      group_by(subclass_by_activity_condition, activity_condition, subclass) |>
      distinct(subclass_by_activity_condition) |> 
      mutate(subclass_by_activity_condition = factor(subclass_by_activity_condition, levels = rownames(plotting_matrix))) |> 
      arrange(subclass_by_activity_condition)
    
    left_annotation <- rowAnnotation(
      activity_condition = df_subclass_annotation$activity_condition,
      subclass = df_subclass_annotation$subclass,
      col = list(
        activity_condition = activity_colors,
        subclass = subclass_colors
      ),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    
    right_annotation <- rowAnnotation(
      subclass = df_subclass_annotation$subclass,
      col = list(
        activity_condition = activity_colors,
        subclass = subclass_colors
      ),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    
    top_annotation <- HeatmapAnnotation(
      classification = df_gene_annotation$classification_TF,
      col = list(
        classification = c(
          # "ERG_0" = "#D81B60", "LRG_0" = "#FFC107", "both_0" = "#82A6B1", 
          "ERG_0" = "gray80", "LRG_0" = "gray80", "both_0" = "gray80", 
          "ERG_1" = "black", "LRG_1" = "black", "both_1" = "black"
          ),
        subclass = subclass_colors
      ),
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
    
    ieg_symbols <- LoadGeneList("IEG")
    bottom_annotation <- HeatmapAnnotation(
      iegs = anno_mark(
        at = which(colnames(plotting_matrix) %in% ieg_symbols),
        labels = intersect(colnames(plotting_matrix), ieg_symbols), 
        side = "bottom"
      ),
      show_annotation_name = FALSE
    )
    
    
    # plot ----------------------------------------------------------------------------
    print(rownames(plotting_matrix))
    p <- Heatmap(
      plotting_matrix, # exclude the first column (subclass_by_activity_condition)
      name = 'HC_average',
      heatmap_legend_param = list(
        title_position = "topcenter",
        direction = "horizontal"
      ),
      # col = circlize::colorRamp2(c(-2, 0, 2), hcl_palette = 'Blue-Red 2'),
      col = circlize::colorRamp2(c(-3,0,3), hcl_palette = 'RdBu', reverse=T),
      row_order = subclass_list,
      cluster_rows = F,
      column_dend_reorder = T,
      top_annotation = top_annotation,
      bottom_annotation = bottom_annotation,
      height = unit(7.5*length(subclass_list), 'mm'),
      show_row_names = F,
      show_column_names = F,
      show_heatmap_legend = F,
      use_raster = F,
    )
    draw(p, heatmap_legend_side = 'bottom')
    plots <- c(plots, list(p)) # add to list of plots
  } # ERG/LRG loop
} # gigaclass loop

# htShiny(plots[[2]], action = 'hover', output_ui_float = T)


# save ----------------------------------------------------------------------------
if (SAVE_PLOTS) {
  print('Saving plots...')
  
  save_dir <- "05-results/NEWT/raw_R_plots/"
  file_str <- "DGE-heatmap_all_subclasses_pseudobulk_"
  file_path <- paste0(save_dir, file_str)  
  # as PNG
  png(paste0(file_path, 'ERG_Ex.png'), width = 8, height = 5, units = "in", res = 900); draw(plots[[1]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_path, 'ERG_In.png'), width = 8, height = 5, units = "in", res = 900); draw(plots[[3]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_path, 'ERG_Gl.png'), width = 8, height = 5, units = "in", res = 900); draw(plots[[5]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_path, 'LRG_Ex.png'), width = 5, height = 5, units = "in", res = 900); draw(plots[[2]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_path, 'LRG_In.png'), width = 5, height = 5, units = "in", res = 900); draw(plots[[4]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_path, 'LRG_Gl.png'), width = 5, height = 5, units = "in", res = 900); draw(plots[[6]], heatmap_legend_side = 'bottom'); dev.off()
  
  # as SVG
  svgsave(plot = plots[[1]], filename = glue('{file_str}_ERG_Ex.svg'), path = save_dir, width = 8, height = 5)
  svgsave(plot = plots[[3]], filename = glue('{file_str}_ERG_In.svg'), path = save_dir, width = 8, height = 5)
  svgsave(plot = plots[[5]], filename = glue('{file_str}_ERG_Gl.svg'), path = save_dir, width = 8, height = 5)
  svgsave(plot = plots[[2]], filename = glue('{file_str}_LRG_Ex.svg'), path = save_dir, width = 5, height = 5)
  svgsave(plot = plots[[4]], filename = glue('{file_str}_LRG_In.svg'), path = save_dir, width = 5, height = 5)
  svgsave(plot = plots[[6]], filename = glue('{file_str}_LRG_Gl.svg'), path = save_dir, width = 5, height = 5)
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))