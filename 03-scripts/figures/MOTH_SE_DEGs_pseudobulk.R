# This script does 2 main things:
#   1. DESeq-based ZT DEG analysis in SE condition
#   2. Plots DEG heatmaps across ZT
# 
# These steps are performed for the 4 largest subclasses:
# CA1, DG, astrocytes, and oligodendrocytes

source("03-scripts/R/seq_functions.R")

library(DESeq2)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(seriation)

USE_SHRUNKEN_FC <- F
PLOT_VOLCANOS <- F


# load data ----------------------------------------
nuclei <- LoadDataset('Dec2024')
seurat_subsets <- list(
  CA1   = subset(nuclei, subclass_name == '016 CA1-ProS Glut' & activity_condition %in% c('SE')),
  DG    = subset(nuclei, subclass_name == '037 DG Glut'       & activity_condition %in% c('SE')),
  Astro = subset(nuclei, subclass_name == '319 Astro-TE NN'   & activity_condition %in% c('SE')),
  Oligo = subset(nuclei, subclass_name == '327 Oligo NN'      & activity_condition %in% c('SE'))
)


# define fxns ----------------------------------------
getUnshrunkenConstrastResults <- function(dds, contrast, threshold) {
  res <- results(dds, contrast = contrast, tidy = T) |> 
    dplyr::rename(gene = row) |> 
    mutate(
      classification = 
        ifelse(log2FoldChange > thresh & padj < 0.05, "upregulated", 
               ifelse(log2FoldChange < -thresh & padj < 0.05, "downregulated", 
                      "no_change")
        )
    ) |> 
    arrange(padj)
  
  print(tibble(res))
  
  rownames(res) <- res$gene
  
  return(tibble(res))
}


getShrunkenConstrastResults <- function(dds, contrast_name, threshold) {
  # get results from Wald test
  res_raw <- results(dds, contrast = contrast_name)

  # filter by shrunken lfc
  res_shrink <- lfcShrink(dds = dds,
                          res = res_raw,
                          contrast = contrast_name,
                          type = 'ashr')|>
    as_tibble(rownames = 'gene') |>
    left_join(select(as_tibble(res_raw, rownames = 'gene'), gene, log2FoldChange, lfcSE),
              by = 'gene',
              suffix = c('.shrink', '.raw'),
              relationship = 'one-to-one') |>
    arrange(padj) |>
    relocate(log2FoldChange.raw, log2FoldChange.shrink, lfcSE.raw, lfcSE.shrink, .after = baseMean) |>
    mutate(
      classification =
        ifelse(  # if upregulated
          log2FoldChange.shrink > threshold & padj < 0.05,
          "upregulated",
          ifelse(  # if downregulated
            log2FoldChange.shrink < -threshold & padj < 0.05,
            "downregulated",
            "no_change" # then no change
          )
        )
    )
  print(res_shrink)
  return(res_shrink)
}


volcanoPlot <- function(res, title, threshold) {
  library(ggrepel)
  library(ggsignif)
  
  point_color <- setNames(c("red", "blue", "black"), 
                          c("upregulated", "downregulated", "no_change"))
  
  p <- ggplot(res) +
    aes(x = log2FoldChange, y = -log10(padj), color = classification) +
    geom_point() +
    geom_text_repel(
      data = subset(res, classification != "no_change"),  # Only label DEGs
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed") +
    scale_color_manual(values = point_color) +
    labs(x = "log2FoldChange", y = "-log10(padj)", title = title) +
    theme(legend.position = "none")
  
  print(p)
  return(p)
}


# run loop ----------------------------------------
plots <- list()
for (subclass in names(seurat_subsets)) {
  nuclei_subclass <- seurat_subsets[[subclass]]
  Idents(nuclei_subclass) <- nuclei_subclass$ZT
  
  
  # prep data ----------------------------------------
  # scramble cell identities for control analysis!!
  # scrambled_indices <- sample(seq_len(nrow(nuclei_subclass@meta.data)))
  # nuclei_subclass@meta.data$sample <- nuclei_subclass@meta.data$sample[scrambled_indices]
  # nuclei_subclass@meta.data$ZT <- nuclei_subclass@meta.data$ZT[scrambled_indices]
  
  # only test genes expressed in 1% of cells
  mat <- nuclei_subclass[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- rownames(mat)
  
  nuclei_subclass$Condition <- factor(nuclei_subclass$activity_condition, levels = c("SE", "EE30m", "EE6h"))
  nuclei_subclass$ZT <- factor(nuclei_subclass$ZT, levels = c("ZT0", "ZT4", "ZT12", "ZT16"))
  nuclei_subclass$sample <- factor(nuclei_subclass$sample)
  
  # for cases where activity_condition == EE6h, set ZT_collection to ZT + 6. Otherwise, use ZT. Please use case_when
  nuclei_subclass@meta.data <- nuclei_subclass@meta.data |> 
    mutate(ZT.collection = case_when(
      (ZT == 'ZT0' & activity_condition == 'EE6h') ~ 'ZT6',
      (ZT == 'ZT4' & activity_condition == 'EE6h') ~ 'ZT10',
      (ZT == 'ZT12' & activity_condition == 'EE6h') ~ 'ZT18',
      (ZT == 'ZT16' & activity_condition == 'EE6h') ~ 'ZT22',
      TRUE ~ ZT
    )) |> 
    mutate(ZT.collection = factor(ZT.collection, 
                                  levels = c('ZT0', 'ZT4', 'ZT6', 'ZT10', 'ZT12', 'ZT16', 'ZT18', 'ZT22')))
  
  pseudobulk_counts <- AggregateExpression(nuclei_subclass, group.by = "sample", features = gene_list, return.seurat = FALSE)$RNA
  
  col_data <- nuclei_subclass@meta.data |> 
    as_tibble() |> 
    distinct(sample, Condition, ZT) |> 
    arrange(ZT, Condition) |> 
    mutate(sample = str_replace(sample, '_', '-')) |> 
    column_to_rownames("sample") |> 
    as.data.frame()
  
  col_data <- col_data[colnames(pseudobulk_counts),]
  col_data$sample <- rownames(col_data)
  col_data
  
  # run DESeq ----------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_counts), 
                                colData = col_data, 
                                design = ~ ZT)
  dds <- DESeq(dds)
  thresh = 0.585
  
  if (USE_SHRUNKEN_FC) {
    res_ZT4_vs_ZT0   <- getShrunkenConstrastResults(dds, c("ZT", "ZT4", "ZT0"), thresh)
    res_ZT12_vs_ZT0  <- getShrunkenConstrastResults(dds, c("ZT", "ZT12", "ZT0"), thresh)
    res_ZT16_vs_ZT0  <- getShrunkenConstrastResults(dds, c("ZT", "ZT16", "ZT0"), thresh)
    res_ZT12_vs_ZT4  <- getShrunkenConstrastResults(dds, c("ZT", "ZT12", "ZT4"), thresh)
    res_ZT16_vs_ZT4  <- getShrunkenConstrastResults(dds, c("ZT", "ZT16", "ZT4"), thresh)
    res_ZT16_vs_ZT12 <- getShrunkenConstrastResults(dds, c("ZT", "ZT16", "ZT12"), thresh)
  } else {
    res_ZT4_vs_ZT0   <- getUnshrunkenConstrastResults(dds, c("ZT", "ZT4", "ZT0"), thresh)
    res_ZT12_vs_ZT0  <- getUnshrunkenConstrastResults(dds, c("ZT", "ZT12", "ZT0"), thresh)
    res_ZT16_vs_ZT0  <- getUnshrunkenConstrastResults(dds, c("ZT", "ZT16", "ZT0"), thresh)
    res_ZT12_vs_ZT4  <- getUnshrunkenConstrastResults(dds, c("ZT", "ZT12", "ZT4"), thresh)
    res_ZT16_vs_ZT4  <- getUnshrunkenConstrastResults(dds, c("ZT", "ZT16", "ZT4"), thresh)
    res_ZT16_vs_ZT12 <- getUnshrunkenConstrastResults(dds, c("ZT", "ZT16", "ZT12"), thresh)
  }
  
  # save results to file
  deg_save_path <- '04-analysis/DEGs/Dec2024_ZT_SE_pseudobulk/'
  
  if (PLOT_VOLCANOS) {
    p <- volcanoPlot(res_ZT4_vs_ZT0, glue("ZT4 vs ZT0: {subclass}"), thresh)
    p <- volcanoPlot(res_ZT12_vs_ZT0, glue("ZT12 vs ZT0: {subclass}"), thresh)
    p <- volcanoPlot(res_ZT16_vs_ZT0, glue("ZT16 vs ZT0: {subclass}"), thresh)
    p <- volcanoPlot(res_ZT12_vs_ZT4, glue("ZT12 vs ZT4: {subclass}"), thresh)
    p <- volcanoPlot(res_ZT16_vs_ZT4, glue("ZT16 vs ZT4: {subclass}"), thresh)
    p <- volcanoPlot(res_ZT16_vs_ZT12, glue("ZT16 vs ZT12: {subclass}"), thresh)
  }
  
  
  # heatmap ----------------------------------------
  # assemble all genes and their contrast of origin
  if (USE_SHRUNKEN_FC) {
    all_genes <- bind_rows(
      res_ZT4_vs_ZT0   |> mutate(contrast = "ZT4 vs ZT0")   |> filter(classification != 'no_change'),
      res_ZT12_vs_ZT0  |> mutate(contrast = "ZT12 vs ZT0")  |> filter(classification != 'no_change'),
      res_ZT16_vs_ZT0  |> mutate(contrast = "ZT16 vs ZT0")  |> filter(classification != 'no_change'),
      res_ZT12_vs_ZT4  |> mutate(contrast = "ZT12 vs ZT4")  |> filter(classification != 'no_change'),
      res_ZT16_vs_ZT4  |> mutate(contrast = "ZT16 vs ZT4")  |> filter(classification != 'no_change'),
      res_ZT16_vs_ZT12 |> mutate(contrast = "ZT16 vs ZT12") |> filter(classification != 'no_change')
    )
    # write_csv(all_genes, glue("{deg_save_path}{subclass}_all_DEGs.csv"))
  } else {
    all_genes <- bind_rows(
      res_ZT4_vs_ZT0   |> mutate(contrast = "ZT4 vs ZT0")   |> filter(classification != 'no_change'),
      res_ZT12_vs_ZT0  |> mutate(contrast = "ZT12 vs ZT0")  |> filter(classification != 'no_change'),
      res_ZT16_vs_ZT0  |> mutate(contrast = "ZT16 vs ZT0")  |> filter(classification != 'no_change'),
      res_ZT12_vs_ZT4  |> mutate(contrast = "ZT12 vs ZT4")  |> filter(classification != 'no_change'),
      res_ZT16_vs_ZT4  |> mutate(contrast = "ZT16 vs ZT4")  |> filter(classification != 'no_change'),
      res_ZT16_vs_ZT12 |> mutate(contrast = "ZT16 vs ZT12") |> filter(classification != 'no_change')
    )
    # write_csv(all_genes, glue("{deg_save_path}{subclass}_all_DEGs.csv"))
  }
  
  # remove duplicates
  gene_list <- all_genes |> 
    distinct(gene) |> 
    pull()
  
  # get DESeq-normalized expression
  counts_mat <- counts(dds, normalized = T)
  gene_mat <- counts_mat[gene_list,] |> 
    t() |> 
    as.data.frame() |> 
    rownames_to_column('sample') |> 
    pivot_longer(cols = -sample, names_to = 'gene', values_to = 'expression') |> 
    mutate(ZT = str_extract(sample, "ZT\\d+")) |> 
    mutate(ZT = factor(ZT, levels = c("ZT0", "ZT4", "ZT12", "ZT16"))) |> 
    group_by(gene, ZT) |> 
    mutate(mean_expression = mean(expression)) |> 
    ungroup() |> 
    group_by(gene) |> 
    mutate(norm_expression = scale(mean_expression)) |>   # Z-SCORE EXPRESSION
    # mutate(norm_expression = mean_expression / max(mean_expression)) |>  # MAX-NORMALIZED EXPRESSION
    select(gene, ZT, norm_expression) |> 
    distinct() |> 
    pivot_wider(names_from = ZT, values_from = norm_expression) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  
  
  # get tau ----------------------------------------
  df_expression_zt <- counts_mat[gene_list,] |> 
    as.data.frame() |> 
    rownames_to_column('gene') |> 
    pivot_longer(cols = -gene, names_to = 'sample', values_to = 'expression') |> 
    separate(sample, into = c('ZT', 'sample')) |> 
    mutate(ZT = factor(ZT, levels = c('ZT0', 'ZT4', 'ZT12', 'ZT16'))) |>
    group_by(gene, ZT) |> 
    dplyr::summarize(expression = mean(expression), .groups = 'drop') |> 
    dplyr::mutate(tau = tau(expression), .by = gene) |> 
    arrange(desc(tau)) |> 
    print()
  
  df_expression_zt |> 
    slice_max(expression, by = gene) |> 
    slice_max(tau, by = ZT) |> 
    arrange(ZT)
  
  
  # seriate matrix ----------------------------------------
  o <- seriate(gene_mat, method = 'Heatmap', seriation_method = 'OLO_average') |> get_order(1)
  gene_mat_ordered <- gene_mat[names(o),]
  hc <- hclust(dist(gene_mat), method = "average")
  dend <- as.dendrogram(hc)
  dend <- dendextend::rotate(dend, order = o)
  
  
  # make tau anno ----------------------------------------
  tau_thresh <- 0.75
  genes_tyssowski <- c(LoadGeneList('tyssowski'), 'Adcy1', 'Adcy8')
  genes_clock <- LoadGeneList('circadian')
  
  anno_colors <- setNames(c('red', 'orange', 'black'), c('high_tau', 'clock', 'tyssowski'))
  anno_colors
  
  # df for tau plot
  df_anno <- df_expression_zt |> 
    group_by(gene, tau) |> 
    summarize() |> 
    mutate(gene = factor(gene, levels = rownames(gene_mat_ordered))) |> 
    arrange(gene) |> 
    mutate(color = ifelse(tau > tau_thresh, anno_colors[1], 
                          ifelse(gene %in% genes_clock, anno_colors[2], 
                                 anno_colors[3]))) |> print()
  
  # get high tau genes
  high_tau_genes <- df_anno |> 
    filter(tau > tau_thresh) |> 
    mutate(text_color = anno_colors[1]) |> 
    select(gene, text_color) |> 
    group_by(gene) |> 
    slice_head() |> 
    mutate(type = 'high_tau') |>
    print()
  
  # add any other genes of interest to list
  clock_genes <- tibble(
    gene = genes_clock,
    text_color = anno_colors[2],
    type = 'clock'
  ) |> 
    filter(!gene %in% high_tau_genes$gene)
  tyssowski_genes <- tibble(
    gene = genes_tyssowski,
    text_color = anno_colors[3],
    type = 'tyssowski'
  ) |> 
    filter(!gene %in% high_tau_genes$gene) |>
    filter(!gene %in% clock_genes$gene)
  
  df_marked_genes <- rbind(high_tau_genes, clock_genes, tyssowski_genes) |> 
    filter(gene %in% rownames(gene_mat_ordered)) |> 
    mutate(gene = factor(gene, levels = rownames(gene_mat_ordered))) |> 
    arrange(gene) |> 
    left_join(read_csv("04-analysis/gene_biotypes.csv", show_col_types = FALSE),
              by = c('gene' = 'gene_name')) |> 
    print(n = 31)
  
  anno_right <- HeatmapAnnotation(
    tau = anno_barplot(
      df_anno$tau, 
      bar_width = 0.7,
      width = unit(1, 'in'),
      ylim = c(0,1),
      gp = gpar(fill = df_anno$color, col = df_anno$color),
      which = 'row'
    ),
    genes = anno_mark(
      at = which(rownames(gene_mat_ordered) %in% df_marked_genes$gene),
      labels = intersect(rownames(gene_mat_ordered), df_marked_genes$gene),
      lines_gp = gpar(col = df_marked_genes$text_color),
      labels_gp = gpar(col = df_marked_genes$text_color, angle = 315),
      side = 'right',
      which = 'row'
    ), 
    which = 'row',
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  
  # find km clusters ----------------------------------------
  # Extract seeds and orders that work well for each subclass.
  # These orders were determined by running kmeans multiple times to
  # find an aesthetically pleasing order.
  subclass_seeds_orders <- tibble(
    subclass_title = c('CA1', 'DG', 'Astro', 'Oligo'),
    seed    = c(17, 1, 9, 5),
    km_order   = list(c(4,1,3,2), c(4,2,3,1), c(1,3,2,4), c(2,4,1,3))
  )
  subclass_specifics <- subclass_seeds_orders |> filter(subclass_title == subclass) |> print()
  seed <- subclass_specifics$seed
  desired_order <- subclass_specifics$km_order[[1]]
  
  # run km
  set.seed(seed)
  km <- kmeans(gene_mat_ordered, centers = 4)
  split <- factor(km$cluster, levels = desired_order)
  
  
  # plot ----------------------------------------
  width_scalar <- 2
  width_params <- list(
    width = unit(width_scalar, 'in'),
    row_dend_width = unit(width_scalar * 0.075, 'in')
    # tau_width = unit(width_scalar * 0., 'in')
  )
  # anno_right@width <- width_params$tau_width
  hm <- Heatmap(
    gene_mat_ordered,
    right_annotation = anno_right,
    row_split = split,
    cluster_rows = T,
    cluster_row_slices = F,  # keeps your factor level order for the slices
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    show_heatmap_legend = F,
    row_title = NULL,
    width = width_params$width, # width of the heatmap body
    row_dend_width = width_params$row_dend_width, # width of the row dendrogram
    col = colorRamp2(c(-2,0,2), hcl_palette = "inferno")
  ) |> draw()
  
  plots <- c(plots, setNames(list(hm), subclass))
}

# InteractiveComplexHeatmap::htShiny(hm)

if (SAVE_PLOTS) {
  for (subclass in names(plots)) {
    subclass_fname <- ShrinkSubclassName(subclass)
    save_path <- "05-results/MOTH/raw_R_plots"
    # png
    png(glue("{save_path}/SE__ZT_heatmap_{subclass_fname}.png"),
        width = 6, height = 10, units = "in", res = 900)
    draw(plots[[subclass]]); dev.off()
    
    # svg
    svgsave(plot = plots[[subclass]],
            path = save_path,
            filename = glue("SE_ZT_heatmap_{subclass_fname}"),
            width = 6, height = 10)
  }
}
