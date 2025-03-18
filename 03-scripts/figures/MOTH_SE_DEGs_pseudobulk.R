# Load required libraries
library(ggplot2)
library(ggrepel)
library(ggsignif)
library(dplyr)
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

source("03-scripts/R/seq_functions.R")


# define fxns -----------------------------------------------------------
getConstrastResults <- function(dds, contrast, threshold) {
  res <- results(dds, contrast = contrast) |> 
    as.data.frame() |> 
    rownames_to_column("gene") |> 
    relocate(gene) |> 
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
subclass_list <- c('016 CA1-ProS Glut', '319 Astro-TE NN', '327 Oligo NN')


for (subclass in subclass_list) {
  celltype <- subset(nuclei, (subclass_name == subclass & activity_condition %in% c('SE')))
  Idents(celltype) <- celltype$ZT
  
  
  # prep data -----------------------------------------------------------
  # scrambled_indices <- sample(seq_len(nrow(celltype@meta.data)))   # scramble cell identities!!
  # celltype@meta.data$sample <- celltype@meta.data$sample[scrambled_indices]
  # celltype@meta.data$ZT <- celltype@meta.data$ZT[scrambled_indices]
  
  # only test genes expressed in 5% of cells
  mat <- celltype[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- rownames(mat)
  
  celltype$Condition <- factor(celltype$activity_condition, levels = c("SE", "EE30m", "EE6h"))
  celltype$ZT <- factor(celltype$ZT, levels = c("ZT0", "ZT4", "ZT12", "ZT16"))
  celltype$sample <- factor(celltype$sample)
  
  # for cases where activity_condition == EE6h, set ZT_collection to ZT + 6. Otherwise, use ZT. Please use case_when
  celltype@meta.data <- celltype@meta.data |> 
    mutate(ZT.collection = case_when(
      (ZT == 'ZT0' & activity_condition == 'EE6h') ~ 'ZT6',
      (ZT == 'ZT4' & activity_condition == 'EE6h') ~ 'ZT10',
      (ZT == 'ZT12' & activity_condition == 'EE6h') ~ 'ZT18',
      (ZT == 'ZT16' & activity_condition == 'EE6h') ~ 'ZT22',
      TRUE ~ ZT
    )) |> 
    mutate(ZT.collection = factor(ZT.collection, 
                                  levels = c('ZT0', 'ZT4', 'ZT6', 'ZT10', 'ZT12', 'ZT16', 'ZT18', 'ZT22')))
  
  pseudobulk_counts <- AggregateExpression(celltype, group.by = "sample", features = gene_list, return.seurat = FALSE)$RNA
  
  col_data <- celltype@meta.data |> 
    as_tibble() |> 
    distinct(sample, Condition, ZT) |> 
    mutate(sample = str_replace(sample, '_', '-')) |> 
    column_to_rownames("sample") |> 
    filter(!is.na(Condition)) |> 
    arrange(ZT, Condition)
  
  col_data <- col_data[colnames(pseudobulk_counts),]
  col_data$sample <- rownames(col_data)
  
  
  # run DESeq -----------------------------------------------------------
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_counts), 
                                colData = col_data, 
                                design = ~ ZT)
  dds <- DESeq(dds)
  
  thresh = 0.585
  
  # ZT4 vs ZT0 ----
  contrast = c("ZT", "ZT4", "ZT0")
  res_ZT4_vs_ZT0 <- getConstrastResults(dds, contrast, thresh)
  p <- volcanoPlot(res_ZT4_vs_ZT0, glue("ZT4 vs ZT0: {subclass}"), thresh)
  
  
  # ZT12 vs ZT0 ----
  contrast = c("ZT", "ZT12", "ZT0")
  res_ZT12_vs_ZT0 <- getConstrastResults(dds, contrast, thresh)
  p <- volcanoPlot(res_ZT12_vs_ZT0, glue("ZT12 vs ZT0: {subclass}"), thresh)
  
  
  # ZT16_vs_ZT0 ----
  contrast = c("ZT", "ZT16", "ZT0")
  res_ZT16_vs_ZT0 <- getConstrastResults(dds, contrast, thresh)
  p <- volcanoPlot(res_ZT16_vs_ZT0, glue("ZT16 vs ZT0: {subclass}"), thresh)
  
  
  # ZT12 vs ZT4 ----
  contrast = c("ZT", "ZT12", "ZT4")
  res_ZT12_vs_ZT4 <- getConstrastResults(dds, contrast, thresh)
  p <- volcanoPlot(res_ZT12_vs_ZT4, glue("ZT12 vs ZT4: {subclass}"), thresh)
  
  
  # ZT16 vs ZT4 ----
  contrast = c("ZT", "ZT16", "ZT4")
  res_ZT16_vs_ZT4 <- getConstrastResults(dds, contrast, thresh)
  volcanoPlot(res_ZT16_vs_ZT4, glue("ZT16 vs ZT4: {subclass}"), thresh)
  
  
  # ZT16 vs ZT12 ----
  contrast = c("ZT", "ZT16", "ZT12")
  res_ZT16_vs_ZT12 <- getConstrastResults(dds, contrast, thresh)
  p <- volcanoPlot(res_ZT16_vs_ZT12, glue("ZT16 vs ZT12: {subclass}"), thresh)
  
  
  # heatmap -----------------------------------------------------------
  # assemble all genes and their contrast of origin
  all_genes <- bind_rows(
    res_ZT4_vs_ZT0 |> mutate(contrast = "ZT4 vs ZT0") |> filter(classification != 'no_change'),
    res_ZT12_vs_ZT0 |> mutate(contrast = "ZT12 vs ZT0") |> filter(classification != 'no_change'),
    res_ZT16_vs_ZT0 |> mutate(contrast = "ZT16 vs ZT0") |> filter(classification != 'no_change'),
    res_ZT12_vs_ZT4 |> mutate(contrast = "ZT12 vs ZT4") |> filter(classification != 'no_change'),
    res_ZT16_vs_ZT4 |> mutate(contrast = "ZT16 vs ZT4") |> filter(classification != 'no_change'),
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
    mutate(z_expression = scale(mean_expression)) |> 
    select(gene, ZT, z_expression) |> 
    distinct() |> 
    pivot_wider(names_from = ZT, values_from = z_expression) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  
  
  # run tau code on ZT genes -----------------------------------------------------------
  df_expression_zt <- AverageExpression(celltype, 
                                        features = gene_list, 
                                        group.by = 'ZT', 
                                        assay = 'SCT', 
                                        layer = 'data')$SCT
  df_zt <- df_expression_zt |> 
    as.data.frame() |> 
    rownames_to_column('gene') |>
    pivot_longer(cols = -gene, names_to = 'zt', values_to = 'expression') |>
    group_by(gene) |> 
    mutate(tau = tau(expression)) |> 
    arrange(desc(tau)) |> 
    print()
  
  df_zt |> 
    group_by(gene) |> 
    slice_max(expression) |> 
    group_by(zt) |> 
    slice_max(tau) |> 
    arrange(desc(tau))
  
  VlnPlot(celltype, features = 'Hif3a')
  
  
  # seriate matrix -----------------------------------------------------------
  o <- seriate(gene_mat, method = 'Heatmap', seriation_method = 'HC_average') |> get_order(1)
  gene_mat_ordered <- gene_mat[names(o),]
  hc <- hclust(dist(gene_mat), method = "average")
  dend <- as.dendrogram(hc)
  dend <- dendextend::rotate(dend, order = o)
  
  # make tau anno -----------------------------------------------------------
  tau_thresh <- 0.75
  df_anno <- df_zt |> 
    slice_head() |> 
    mutate(gene = factor(gene, levels = rownames(gene_mat))) |> 
    mutate(color = ifelse(tau > tau_thresh, 'red', 'black')) |> 
    arrange(gene)
  
  # get high tau genes
  high_tau_genes <- df_anno |> 
    filter(tau > tau_thresh) |> 
    mutate(text_color = 'red') |> 
    select(gene, text_color)
  
  # add any other genes of interest to list
  curated_genes <- tibble(
    gene = c('Per1', 'Per2', 'Bmal1', 'Cry1', 'Cry2', 'Adcy1', 'Adcy8'),
    text_color = 'black'
  )
  
  df_marked_genes <- rbind(high_tau_genes, curated_genes) |> 
    mutate(gene = factor(gene, levels = rownames(gene_mat))) |> 
    arrange(gene)
  
  anno_right <- HeatmapAnnotation(
    tau = anno_barplot(
      df_anno$tau, 
      bar_width = 0.7,
      width = unit(3, 'cm'),
      ylim = c(0,1),
      gp = gpar(fill = df_anno$color, col = df_anno$color),
      which = 'row'
    ),
    genes = anno_mark(
      at = which(rownames(gene_mat) %in% df_marked_genes$gene),
      labels = intersect(rownames(gene_mat), df_marked_genes$gene),
      lines_gp = gpar(col = df_marked_genes$text_color),
      labels_gp = gpar(col = df_marked_genes$text_color, angle = 315),
      side = 'right',
      which = 'row'
    ), 
    which = 'row',
    show_annotation_name = FALSE,
    show_legend = FALSE
  )
  
  
  # plot -----------------------------------------------------------
  hm <- Heatmap(
    gene_mat,
    col = colorRamp2(c(-2,0,2), hcl_palette = "blue-red2"),
    column_names_gp = gpar(rot = 180),
    right_annotation = anno_right,
    cluster_columns = FALSE,
    cluster_rows = dend,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_heatmap_legend = FALSE
  )
  
  if (SAVE_PLOTS) {
    subclass_fname <- ShrinkSubclassName(subclass)
    png(glue("05-results/figureZT/raw_R_plots/SE_ZT_heatmap_{subclass_fname}.png"), 
        width = 5, height = 10, units = "in", res = 900)
    draw(hm)
    dev.off()
  } else{
    draw(hm)
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