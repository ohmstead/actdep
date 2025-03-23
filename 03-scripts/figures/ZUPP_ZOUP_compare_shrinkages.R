# This script will use DESeq2 to run an interaction analysis on whether genes
# are signficiantly changing in expression at different times of day.

library(ggplot2)
library(dplyr)
library(glue)

library(Seurat)
library(DESeq2)
library(ComplexHeatmap)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset("Dec2024")
zt_colors <- LoadZTColors()
activity_condition <- LoadActivityColors()

# establish subclasses to use ------------------------------------
subclass_list <- LoadSubclassesToUse(nuclei)
subclasses <- subclass_list[c(1,4, 17, 19)]

seurat_subsets <- list(
  ca1   = nuclei |> subset(subclass_name == subclass_list[1]  & activity_condition %in% c('SE', 'EE30m', 'EE6h')),
  dg    = nuclei |> subset(subclass_name == subclass_list[4]  & activity_condition %in% c('SE', 'EE30m', 'EE6h')),
  astro = nuclei |> subset(subclass_name == subclass_list[17] & activity_condition %in% c('SE', 'EE30m', 'EE6h')),
  oligo = nuclei |> subset(subclass_name == subclass_list[19] & activity_condition %in% c('SE', 'EE30m', 'EE6h'))
)


for (subclass in names(seurat_subsets)) {
# make DESeq -----------------------------------------------------
  # only test genes expressed in 5% of cells
  nuclei_subclass <- seurat_subsets[[subclass]]
  mat <- nuclei_subclass[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- rownames(mat)
  
  nuclei_subclass$activity_condition <- factor(nuclei_subclass$activity_condition, levels = c("SE", "EE30m", "EE6h"))
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
  pseudobulk_counts <- AggregateExpression(nuclei_subclass, 
                                           group.by = "sample", 
                                           features = gene_list, 
                                           return.seurat = FALSE)$RNA
  col_data <- nuclei_subclass@meta.data |> 
    as_tibble() |> 
    distinct(sample, activity_condition, ZT) |> 
    mutate(sample = str_replace(sample, '_', '-')) |> 
    column_to_rownames("sample") |> 
    filter(!is.na(activity_condition)) |> 
    arrange(ZT, activity_condition)
  
  col_data <- col_data[colnames(pseudobulk_counts),]
  col_data$sample <- rownames(col_data)
  
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_counts), 
                                colData = col_data, 
                                design = ~ activity_condition + ZT + activity_condition:ZT)
  
  
  # run DESeq --------------------------------------------
  dds <- DESeq(dds)
  FC_thresh <- 0.585
  
  # look at results
  con <- c("activity_condition", "EE30m", "SE")
  res <- results(dds, contrast = con, alpha = 0.05)
  res_shrink <- lfcShrink(dds, coef = 'activity_condition_EE30m_vs_SE', res=res)
  
  tmp_shrink <- res_shrink |> 
    as.data.frame() |> 
    rownames_to_column('gene') |> 
    as_tibble() |> 
    # filter(padj < 0.05, abs(log2FoldChange) > 0.585) |>
    arrange(desc(log2FoldChange)) |>
    print()
  
  tmp <- res |> 
    as.data.frame() |>
    rownames_to_column('gene') |> 
    as_tibble() |> 
    # filter(padj < 0.05, abs(log2FoldChange) > FC_thresh) |>
    arrange(desc(log2FoldChange)) |> 
    print()
  
  setdiff(tmp$gene, tmp_shrink$gene)
  
  # make df to compare shrunken and unshrunken log2FC
  df <- tmp |> 
    left_join(tmp_shrink, by = 'gene', suffix = c('.unshrink', '.shrink')) |>
    print()
  
  
  df <- df |> 
    filter(baseMean.unshrink < 500) |>
    mutate(color = ifelse(padj.unshrink < 0.05, 'red', 'black')) |>
    mutate(color = ifelse(gene %in% LoadGeneList(), 'green', color))
  
  ggplot(df) +
    aes(x = baseMean.unshrink, 
        y = log2FoldChange.unshrink - log2FoldChange.shrink, 
        color = log2FoldChange.unshrink) +
    geom_point(alpha = 1) +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
}