# load libs and data ----
print("Loading libraries and data...")
library(Seurat)
library(tidyverse)
library(readxl)
library(glue)
library(DESeq2)
library(patchwork)
library(seriation)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)

source('03-scripts/R/seq_functions.R')
activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# establish subclasses to use ----------------------------------------------------
print('Subsetting Seurat object...')
subclass_list <- LoadSubclassesToUse(nuclei, ascertainment = 'custom')
subclass_sets <- LoadSubclassesToUse(nuclei, as_gigaclasses = T)
seurat_subsets <- LoadDataset('Dec2024', as_gigaclass = T)


# loop thru sets of subclasses
plots <- list()
list_df_expression <- list(excitatory = NULL, inhibitory = NULL, glia = NULL)
for (gigaclass in names(subclass_sets)) {
  subclass_list <- subclass_sets[[gigaclass]]
  
  # get list of subclasses to use
  nuclei_subclass <- seurat_subsets[[gigaclass]]
  Idents(nuclei_subclass) <- nuclei_subclass$subclass_name
  
  
  # read in DEGs ----------------------------------------------------
  print("Getting DEGs for all subclasses and contrasts...")
  contrast_list <- c("EE30m_vs_SE", "EE6h_vs_SE")
  dir_deg <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
  deg_files <- list.files(dir_deg, full.names = TRUE)
  
  # loop thru CSVs and load in DEGs
  df_all_DEGs <- data.frame()  # init empty df
  for (file in deg_files) {
    fname <- basename(file)
    fname <- str_remove(fname, "\\.csv$")
    # Extract the contrast (last two underscore-separated parts in filename)
    contrast <- str_extract(fname, "[^_]+_[^_]+_[^_]+$")
    # Extract the subclass (everything before the contrast in filename)
    subclass <- str_remove(fname, paste0("_", contrast, "$"))
    subclass <- str_replace_all(subclass, "_", " ")
    subclass <- str_sub(subclass, 1, -2)
    
    # skip excluded subclasses or contrasts
    if (!(contrast %in% contrast_list)) {
      next
    } else if (!(subclass %in% subclass_list)) {
      next
    }
    
    # add gene list to df
    deg <- read_csv(file) |> 
      arrange(desc(log2FoldChange)) |> 
      filter(classification != 'no_change')
      # filter(padj < 0.05)
    
    deg$subclass <- subclass
    deg$contrast <- contrast
    
    df_all_DEGs <- rbind(df_all_DEGs, deg)
  }
  
  
  # classify DEGs ----------------------------------------------------------------------------
  print("Filtering DEGs")
  
  df_gene_classifications <- df_all_DEGs |> 
    mutate(is_ERG = contrast == "EE30m_vs_SE",
           is_LRG = contrast == 'EE6h_vs_SE') |> 
    group_by(gene) |>
    summarize(
      is_ERG = any(is_ERG, na.rm = T),
      is_LRG = any(is_LRG, na.rm = T),
    ) |> 
    mutate(classification = case_when(
      is_ERG & is_LRG ~ "both",
      is_ERG ~ "ERG",
      is_LRG ~ "LRG",
      TRUE ~ NA_character_
    )) |> print()
  
  # for each subclass, get expression for each activity_condition
  df_distinct_DEGs <- df_all_DEGs |> 
    filter(subclass %in% subclass_list) |> 
    filter(abs(log2FoldChange) >= 0.585) |> 
    mutate(activity_condition = case_when(
      contrast == "EE30m_vs_SE" ~ "EE30m",
      contrast == "EE6h_vs_SE" ~ "EE6h",
      contrast == "KA30m_vs_SE" ~ "KA30m",
      contrast == "KA6h_vs_SE" ~ "KA6h",
      TRUE ~ NA_character_
    )) |>
    filter(!str_detect(contrast, "KA")) |>
    group_by(gene, subclass, contrast) |> 
    summarize(log2FoldChange = log2FoldChange, .groups = 'drop') |> 
    left_join(df_gene_classifications, by = 'gene') |> 
    mutate(subclass_by_contrast = paste(subclass, contrast, sep = " x ")) |> 
    distinct(gene, .keep_all = TRUE) |> 
    mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) |> 
    select(-subclass) # remove subclass for future joins
  
  
  # normalize data ----------------------------------------------------------------------------
  print("Normalizing expression data...")
  # DESeq log2FC
  pseudobulk_counts <- AggregateExpression(
      nuclei_subclass,
      assays = 'RNA',
      features = df_distinct_DEGs$gene,
      group.by = c("subclass_name", "activity_condition")
    )$RNA |>
    as.matrix()
  
  coldata <- data.frame(
    group = colnames(pseudobulk_counts),
    stringsAsFactors = FALSE
  )
  coldata <- coldata |>
    separate(group, into = c("subclass", "activity_condition"), sep = "_", remove = F)
  rownames(coldata) <- coldata$group
  
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_counts,
    colData = coldata,
    design = ~ subclass + activity_condition
  )
  
  # Set "SE" as the reference level for activity_condition:
  dds$activity_condition <- relevel(dds$activity_condition, ref = "SE")
  dds <- DESeq(dds)
  
  # log-transform normalized counts
  norm_counts <- counts(dds, normalized = TRUE)
  norm_log2 <- log2(norm_counts + 1)  # Adding a pseudocount to avoid log(0)
  df_norm <- as.data.frame(norm_log2) |>
    rownames_to_column(var = "gene") |>
    pivot_longer(cols = -gene, names_to = "group", values_to = "avg_expression_log2")
  
  # re-extract variables
  df_norm <- df_norm |>
    separate(group, into = c("subclass", "activity_condition"), sep = "_") |>
    mutate(subclass = str_sub(subclass, 2, -1))
  
  # calculate log2FoldChange for each gene x subclass combo
  df_expression.DESeq <- df_norm |>
    filter(!str_detect(activity_condition, 'KA')) |>
    group_by(gene, subclass) |>
    mutate(log2FC_calculated = avg_expression_log2 - avg_expression_log2[activity_condition == "SE"]) |>
    left_join(df_distinct_DEGs, by = 'gene') |>
    group_by(gene, subclass) |>
    mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
    mutate(activity_condition = factor(activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))) |>
    mutate(subclass = factor(subclass, levels = subclass_list)) |>
    relocate(gene, subclass, activity_condition, log2FC_calculated) |>
    ungroup()

  
  # for genes in the "both" category, force ERG/LRG ----------------------------------------------------------------------------
  df_expression <- df_expression.DESeq
  expression_range <- df_expression |>
    group_by(gene, subclass) |>
    summarise(
      max_expr_condition = activity_condition[which.max(log2FC_calculated)],
      min_expr_condition = activity_condition[which.min(log2FC_calculated)]
    )
  
  # re-classify as ERG or LRG
  df_expression <- df_expression |>
  left_join(expression_range, by = c("gene", "subclass")) |>
  mutate(classification = if_else(
    classification == "both",
    if_else(
      direction == "up",
      if_else(max_expr_condition == "EE30m", "ERG", "LRG"), # if upregulated
      if_else(min_expr_condition == "EE30m", "ERG", "LRG"),   # if downregulated
      ),
    classification  # For genes not labeled "both", keep original classification.
    ),
  )
  
  # for any genes that may be ERG in one subclass and LRG in another, force to modal classification
  df_expression <- df_expression |>
    group_by(gene) |>
    mutate(modal_classification = names(sort(table(classification), decreasing = TRUE))[1]) |>
    ungroup() |>
    mutate(classification = modal_classification) |>
    select(-modal_classification)
  
  # store 
  list_df_expression[[gigaclass]] <- df_expression
  
  # get TF genes ----------------------------------------------------------------------------
  print("Make mat for heatmap...")
  
  # get gene classifications from df_expression
  df_gene_levels <- df_distinct_DEGs |>
    select(-classification) |> # remove old classification that includes 'both' category
    left_join(distinct(select(df_expression, gene, classification))) |> 
    filter(!is.na(classification))
  
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
    select(gene, subclass_by_activity_condition, log2FC_calculated) |> 
    pivot_wider(names_from = 'gene', values_from = 'log2FC_calculated') |> 
    column_to_rownames('subclass_by_activity_condition') |> 
    as.matrix()
  
  
  # seriate gene columns ----------------------------------------------------------------------------
  print('Seriating gene order')
  expression_matrix <- expression_matrix[,df_gene_levels$gene]
  
  erg_matrix <- expression_matrix[str_detect(rownames(expression_matrix), 'EE30m'), df_gene_levels$classification == 'ERG']
  lrg_matrix <- expression_matrix[str_detect(rownames(expression_matrix), 'EE6h'),  df_gene_levels$classification == 'LRG']
  
  ser_method <- 'OLO_average'
  
  erg_seriation <- seriate(erg_matrix, method = 'Heatmap', seriation_method = ser_method)
  lrg_seriation <- seriate(lrg_matrix, method = 'Heatmap', seriation_method = ser_method)
  erg_order <- erg_seriation |> get_order(2)
  lrg_order <- lrg_seriation |> get_order(2)
  col_order <- c(erg_order, lrg_order+length(erg_order))
  
  erg_matrix <- erg_matrix[,names(erg_order)]
  lrg_matrix <- lrg_matrix[,names(lrg_order)]
  
  mat <- list(erg_matrix, lrg_matrix)
  
  # relevel/arrange df_gene_levels for annotation objects
  df_gene_levels <- df_gene_levels |> 
    mutate(gene = factor(gene, levels = names(col_order))) |> 
    arrange(gene)
  
  
  # make annotation objects ----------------------------------------------------------------------------
  print('Making annotation objects and plotting')
  df_looping <- data.frame(
    gene_classification = c('ERG', 'LRG'),
    condition = c('EE30m', 'EE6h'),
    seriation = I(list(erg_seriation, lrg_seriation))
    )
  for (k in 1:2) { # loop ERG/LRG figures
    looping <- df_looping[k,]
    plotting_matrix <- mat[[k]]
    
    df_gene_annotation <- df_gene_levels |> 
      filter(classification == looping$gene_classification) |> 
      mutate(gene = factor(gene), levels = colnames(plotting_matrix)) |> 
      arrange(gene)
      
    df_subclass_annotation <- df_expression |>
      filter(classification == looping$gene_classification) |>
      filter(activity_condition == looping$condition) |> 
      group_by(subclass_by_activity_condition, activity_condition, subclass) |>
      distinct(subclass_by_activity_condition) |> 
      mutate(subclass_by_activity_condition = factor(subclass_by_activity_condition, levels = rownames(expression_matrix))) |> 
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
          "ERG_0" = "#D81B60", "LRG_0" = "#FFC107", "both_0" = "#82A6B1", 
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
        at = which(df_gene_annotation$gene %in% ieg_symbols),
        labels = intersect(colnames(plotting_matrix), ieg_symbols), 
        side = "bottom"
      ),
      show_annotation_name = FALSE
    )
    
    
    # plot ----------------------------------------------------------------------------
    p <- Heatmap(
      plotting_matrix, # exclude the first column (subclass_by_activity_condition)
      name = 'HC_average',
      heatmap_legend_param = list(
        title_position = "topcenter",
        direction = "horizontal"
      ),
      col = circlize::colorRamp2(c(-2, 0, 2), hcl_palette = 'Blue-Red 2'),
      cluster_rows = FALSE,
      # cluster_columns = as.dendrogram(looping$seriation[[1]][[2]]),
      column_dend_reorder = TRUE,
      # cluster_columns = F,
      # left_annotation = left_annotation,
      # right_annotation = right_annotation,
      top_annotation = top_annotation,
      bottom_annotation = bottom_annotation,
      row_gap = unit(2, "mm"),
      column_gap = unit(2, "mm"),
      show_row_names = FALSE,
      show_column_names = TRUE,
      use_raster = FALSE
    )
    draw(p, heatmap_legend_side = 'bottom')
    plots <- c(plots, list(p)) # add to list of plots
  } # end loop for ERG/LRG figures
} # end loop for this set of subclasses
htShiny(plots[[1]], action = 'hover', output_ui_float = T)


# save ----------------------------------------------------------------------------
if (SAVE_PLOTS) {
  file_str <- "05-results/figure1/raw_R_plots/DGE-heatmap_all_subclasses_pseudobulk_"
  #ERG plots
  png(paste0(file_str, '1.png'), width = 8, height = 5, units = "in", res = 900); draw(plots[[1]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_str, '3.png'), width = 8, height = 5, units = "in", res = 900); draw(plots[[3]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_str, '5.png'), width = 8, height = 5, units = "in", res = 900); draw(plots[[5]], heatmap_legend_side = 'bottom'); dev.off()
  # LRG plots
  png(paste0(file_str, '2.png'), width = 5, height = 5, units = "in", res = 900); draw(plots[[2]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_str, '4.png'), width = 5, height = 5, units = "in", res = 900); draw(plots[[4]], heatmap_legend_side = 'bottom'); dev.off()
  png(paste0(file_str, '6.png'), width = 5, height = 5, units = "in", res = 900); draw(plots[[6]], heatmap_legend_side = 'bottom'); dev.off()
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))