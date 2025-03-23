# Load required libraries
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(tidyr)
library(tibble)
library(stringr)
library(glue)

library(Seurat)
library(DESeq2)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(seriation)
library(dplyr)

source("03-scripts/R/seq_functions.R")


# define fxns -----------------------------------------------------------
getConstrastResults <- function(dds, contrast, threshold) {
  res <- results(dds, contrast = contrast, tidy = T) |> 
    dplyr::rename(gene = row) |> 
    mutate(
      classification = 
        ifelse(log2FoldChange > thresh & padj < 0.05, "upregulated", 
               ifelse(log2FoldChange < -thresh & padj < 0.05, "downregulated", 
                      "no_change")
        )
    )
  
  print(tibble(res))
  
  rownames(res) <- res$gene
  
  return(tibble(res))
}


volcanoPlot <- function(res, title, threshold) {
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


data_summary <- function(data, varname, groupnames) {
  #+++++++++++++++++++++++++
  # Function to calculate the mean and the standard deviation
  # for each group
  #+++++++++++++++++++++++++
  # data : a data frame
  # varname : the name of a column containing the variable
  #to be summariezed
  # groupnames : vector of column names to be used as
  # grouping variables
  require(plyr)
  
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE) / sqrt(length(x[[col]])) )
  }
  
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# load data ----------------------------------------------------------------------------------------------------
nuclei <- LoadDataset('Dec2024')
seurat_subsets <- list(
  CA1   = subset(nuclei, subclass_name == '016 CA1-ProS Glut' & activity_condition %in% c('SE')),
  DG    = subset(nuclei, subclass_name == '037 DG Glut'       & activity_condition %in% c('SE')),
  Astro = subset(nuclei, subclass_name == '319 Astro-TE NN'   & activity_condition %in% c('SE')),
  Oligo = subset(nuclei, subclass_name == '327 Oligo NN'      & activity_condition %in% c('SE'))
)


for (subclass in names(seurat_subsets)) {
  nuclei_subclass <- seurat_subsets[[subclass]]
  Idents(nuclei_subclass) <- nuclei_subclass$ZT
  
  
  # prep data -----------------------------------------------------------
  # scrambled_indices <- sample(seq_len(nrow(nuclei_subclass@meta.data)))   # scramble cell identities!!
  # nuclei_subclass@meta.data$sample <- nuclei_subclass@meta.data$sample[scrambled_indices]
  # nuclei_subclass@meta.data$ZT <- nuclei_subclass@meta.data$ZT[scrambled_indices]
  
  # only test genes expressed in 5% of cells
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
  
  # run DESeq -----------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_counts), 
                                colData = col_data, 
                                design = ~ ZT)
  dds <- DESeq(dds)
  thresh = 0.585
  
  res_ZT4_vs_ZT0 <- getConstrastResults(dds, c("ZT", "ZT4", "ZT0"), thresh)
  res_ZT12_vs_ZT0 <- getConstrastResults(dds, c("ZT", "ZT12", "ZT0"), thresh)
  res_ZT16_vs_ZT0 <- getConstrastResults(dds, c("ZT", "ZT16", "ZT0"), thresh)
  res_ZT12_vs_ZT4 <- getConstrastResults(dds, c("ZT", "ZT12", "ZT4"), thresh)
  res_ZT16_vs_ZT4 <- getConstrastResults(dds, c("ZT", "ZT16", "ZT4"), thresh)
  res_ZT16_vs_ZT12 <- getConstrastResults(dds, c("ZT", "ZT16", "ZT12"), thresh)
  
  # p <- volcanoPlot(res_ZT4_vs_ZT0, glue("ZT4 vs ZT0: {subclass}"), thresh)
  # p <- volcanoPlot(res_ZT12_vs_ZT0, glue("ZT12 vs ZT0: {subclass}"), thresh)
  # p <- volcanoPlot(res_ZT16_vs_ZT0, glue("ZT16 vs ZT0: {subclass}"), thresh)
  # p <- volcanoPlot(res_ZT12_vs_ZT4, glue("ZT12 vs ZT4: {subclass}"), thresh)
  # p <- volcanoPlot(res_ZT16_vs_ZT4, glue("ZT16 vs ZT4: {subclass}"), thresh)
  # p <- volcanoPlot(res_ZT16_vs_ZT12, glue("ZT16 vs ZT12: {subclass}"), thresh)
  
  
  # heatmap -----------------------------------------------------------
  # assemble all genes and their contrast of origin
  all_genes <- bind_rows(
    res_ZT4_vs_ZT0 |> mutate(contrast = "ZT4 vs ZT0")     |> filter(classification != 'no_change'),
    res_ZT12_vs_ZT0 |> mutate(contrast = "ZT12 vs ZT0")   |> filter(classification != 'no_change'),
    res_ZT16_vs_ZT0 |> mutate(contrast = "ZT16 vs ZT0")   |> filter(classification != 'no_change'),
    res_ZT12_vs_ZT4 |> mutate(contrast = "ZT12 vs ZT4")   |> filter(classification != 'no_change'),
    res_ZT16_vs_ZT4 |> mutate(contrast = "ZT16 vs ZT4")   |> filter(classification != 'no_change'),
    res_ZT16_vs_ZT12 |> mutate(contrast = "ZT16 vs ZT12") |> filter(classification != 'no_change')
  )
  
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
    mutate(z_expression = scale(mean_expression)) |>   # Z-SCORE EXPRESSION
    select(gene, ZT, z_expression) |> 
    distinct() |> 
    pivot_wider(names_from = ZT, values_from = z_expression) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  
  
  # get tau -----------------------------------------------------------
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
  
  
  # seriate matrix -----------------------------------------------------------
  o <- seriate(gene_mat, method = 'Heatmap', seriation_method = 'OLO_average') |> get_order(1)
  gene_mat_ordered <- gene_mat[names(o),]
  hc <- hclust(dist(gene_mat), method = "average")
  dend <- as.dendrogram(hc)
  dend <- dendextend::rotate(dend, order = o)
  
  # make tau anno -----------------------------------------------------------
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
  # extract seeds and orders that work well for each subclass
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
  hm <- Heatmap(
    gene_mat_ordered,
    right_annotation = anno_right,
    row_split = split,
    cluster_rows = T,
    cluster_row_slices = F,  # keeps your factor level order for the slices
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    show_heatmap_legend = T,
    row_title = NULL,
    width = unit(4, 'in'), # width of the heatmap body
    row_dend_width = unit(0.3, 'in'),
    col = colorRamp2(c(-2,0,2), hcl_palette = "blue-red2")
  ) |> draw()
  
  if (SAVE_PLOTS) {
    subclass_fname <- ShrinkSubclassName(subclass)
    save_path <- "05-results/MOTH/raw_R_plots"
    # png
    png(glue("{save_path}/SE_ZT_heatmap_{subclass_fname}.png"),
        width = 8, height = 10, units = "in", res = 900)
    draw(hm); dev.off()
    
    # svg
    # svg(glue("{save_path}/SE_ZT_heatmap_{subclass_fname}.svg"),
    svg(glue("{save_path}/SE_ZT_heatmap_LEGEND.svg"),
        width = 8, height = 10, bg = 'transparent')
    draw(hm)
    dev.off()
  }
}

InteractiveComplexHeatmap::htShiny(hm)

# use compareCluster from clusterProfiler ----
# make ensemble biomart object
# ensembl <- useEnsembl(biomart = "ensembl", 
#                       dataset = "mmusculus_gene_ensembl",
#                       version = 113)
# 
# df_background_genes <- getBM(
#   attributes = c("mgi_symbol", "entrezgene_id", "ensembl_gene_id"),
#   mart = ensembl
# )
# df_background_genes$entrezgene_id <- as.character(df_background_genes$entrezgene_id)
# 
# list_cluster_genes <- list(
#   cluster.1 = df_background_genes[df_background_genes$mgi_symbol %in% cluster_genes[[1]], "entrezgene_id"],
#   cluster.2 = df_background_genes[df_background_genes$mgi_symbol %in% cluster_genes[[2]], "entrezgene_id"],
#   cluster.3 = df_background_genes[df_background_genes$mgi_symbol %in% cluster_genes[[3]], "entrezgene_id"],
#   cluster.4 = df_background_genes[df_background_genes$mgi_symbol %in% cluster_genes[[4]], "entrezgene_id"]
# )
# list_cluster_genes <- lapply(list_cluster_genes, as.character)
# 
# ego <- compareCluster(list_cluster_genes, 
#                fun = "enrichGO",
#                ont = 'BP',
#                pvalueCutoff = 0.05,
#                OrgDb = org.Mm.eg.db,
#                universe = df_background_genes$entrezgene_id,
#                )
# ego <- setReadable(ego, org.Mm.eg.db, keyType = 'ENTREZID') |> 
#   as_tibble()
# 
# ego |> 
#   group_by(Cluster) |> 
#   arrange(desc(FoldEnrichment)) |> 
#   head()
# 
# 
# # GO analysis ----
# cluster_genes <- lapply(row_order(hm), function(idx) rownames(mat)[idx])
# background_genes <- rownames(mat)
# 
# 
# terms_1 <- RunGOEnrichment(cluster_genes[[1]]) |> mutate(cluster = 1) |> arrange(p.adjust)
# terms_2 <- RunGOEnrichment(cluster_genes[[2]]) |> mutate(cluster = 2) |> arrange(p.adjust)
# terms_3 <- RunGOEnrichment(cluster_genes[[3]]) |> mutate(cluster = 3) |> arrange(p.adjust)
# terms_4 <- RunGOEnrichment(cluster_genes[[4]]) |> mutate(cluster = 4) |> arrange(p.adjust)
# 
# # remove redundant IDs from each set of top terms
# terms_1_curated <- terms_1 |> slice(-c(4, 5, 6))
# terms_2_curated <- terms_2 |> slice(-c(3, 4))
# terms_3_curated <- terms_3 |> slice(-c(2, 3, 4))
# terms_4_curated <- terms_4 |> slice(-c(5))
# 
# 
# go_terms_all <- rbind(
#   terms_1_curated,
#   terms_2_curated,
#   terms_3_curated,
#   terms_4_curated
# )
# 
# # get the top 5 terms from each cluster
# top5_each_cluster <- go_terms_all |> 
#   arrange(p.adjust) |> 
#   slice_head(n = 5, by = cluster) |> 
#   arrange(cluster)
# 
# # get values in all clusters from top5_each_cluster
# cluster1_values <- terms_1_curated |> right_join(top5_each_cluster, by = "ID", suffix = c('.unique', '.ascertainment'))
# cluster2_values <- terms_2_curated |> right_join(top5_each_cluster, by = "ID", suffix = c('.unique', '.ascertainment'))
# cluster3_values <- terms_3_curated |> right_join(top5_each_cluster, by = "ID", suffix = c('.unique', '.ascertainment'))
# cluster4_values <- terms_4_curated |> right_join(top5_each_cluster, by = "ID", suffix = c('.unique', '.ascertainment'))
# 
# all_cluster_vals <- rbind(
#   cluster1_values, 
#   cluster2_values, 
#   cluster3_values, 
#   cluster4_values
# ) |> 
#   mutate(cluster.ascertainment = factor(cluster.ascertainment)) |> 
#   mutate(Description.unique = factor(Description.unique, levels = unique(top5_each_cluster$Description))) |> 
#   mutate(Description.ascertainment = factor(Description.ascertainment, levels = unique(top5_each_cluster$Description)))
# 
# # plot with geom_tile
# ggplot(all_cluster_vals) +
#   aes(x = cluster.unique, y = Description.ascertainment, fill = -log10(p.adjust.unique)) +
#   geom_tile(color = 'black', linewidth = 0.25) +
#   scale_fill_gradient(limits = c(0,3), oob = scales::squish) +
#   scale_y_discrete(position = 'right', limits = rev) +
#   theme_minimal() +
#   labs(fill = 'z-score') +
#   theme(
#     axis.title.x = element_blank(),
#     axis.title.y = element_blank(),
#     axis.text.x = element_blank(),
#     axis.text.y.right = element_text(size = 12),
#     legend.title = element_blank(),
#     legend.position = 'bottom',
#     legend.title.position = 'top'
#     )
# 