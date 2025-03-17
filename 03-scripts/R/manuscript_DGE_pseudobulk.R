# load libs & data ----
library(Seurat)
library(tidyverse)
library(DESeq2)

source("03-scripts/R/seq_functions.R")

nuclei <- LoadDataset("Dec2024")

csv_write_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"


# findDEG() ----
findDEG <- function(seurat_obj, subclass, group1, group2, logFC_threshold=0.25, pct = 0.01) {
  # print report statement
  print(glue("{subclass} between {group1} and {group2}."))

  # skip is number of cells < 30 per group
  if (sum(seurat_obj$activity_condition == group1) < 30 | 
      sum(seurat_obj$activity_condition == group2) < 30) {
    print(glue("Skipping {subclass} between {group1} and {group2} due to low cell count."))

    return(data.frame())
  }
  
  # extract the contrast
  
  return(DEGs)
}

getConstrastResults <- function(dds, contrast, threshold) {
  res <- results(dds, contrast = contrast) |> 
    as.data.frame() |> 
    rownames_to_column("gene") |> 
    relocate(gene) |> 
    mutate(
      classification = 
        ifelse(log2FoldChange > threshold & padj < 0.05, "upregulated", 
               ifelse(log2FoldChange < -threshold & padj < 0.05, "downregulated", 
                      "no_change")
        )
    ) |> 
    arrange(padj)
    
    print(tibble(res))

  rownames(res) <- res$gene
  
  return(tibble(res))
}


# run DEG analysis for each subclass ----
subclass_list <- LoadSubclassesToUse(nuclei)

log2FC_threshold <- 0.585  # corresponding to 50% change

for (i in seq_along(subclass_list)) {
  subclass <- subclass_list[i]
  print(glue("subclass {i}/{length(subclass_list)} -- {subclass}"))

  # subset nuclei for current subclass
  nuclei_subclass <- subset(nuclei, subclass_name == subclass)
  Idents(nuclei_subclass) <- nuclei_subclass$activity_condition
  subclass_fname <- ShrinkSubclassName(subclass)
  
  # run DESeq for this subclass ----
  # only test genes expressed in 1% of cells
  mat <- nuclei_subclass[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- rownames(mat)
  
  # use DESeq to find DEGs
  pseudobulk_counts <- AggregateExpression(nuclei_subclass, 
                                           group.by = 'sample', 
                                           features = gene_list,
                                           return.seurat = FALSE)$RNA
  
  # put together column data for DESeq  
  col_data <- nuclei_subclass@meta.data |> 
    as_tibble() |> 
    distinct(sample, activity_condition) |> 
    mutate(sample = str_replace(sample, '_', '-')) |> 
    column_to_rownames("sample") |> 
    arrange(activity_condition)
  
  col_data$sample <- rownames(col_data)
  # reorder col_data rows to match pseudobulk_counts
  col_data <- col_data[colnames(pseudobulk_counts), ]
  
  
  # run DESeq ----
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_counts), 
                                colData = col_data, 
                                design = ~activity_condition)
  dds <- DESeq(dds)
  
  # get contrasts available for this subclass
  df_contrasts <- GetSubclassContrasts(nuclei, subclass, cell_cutoff = 30)
  
  # run contrasts for cell groups that meet cell_cutoff
  for (k in seq_along(df_contrasts$group1)) {
    group1 <- as.character(df_contrasts$group1[k])
    group2 <- as.character(df_contrasts$group2[k])
    
    specific_contrast <- c("activity_condition", group1, group2)

    DGE_results <- getConstrastResults(dds, contrast = specific_contrast, log2FC_threshold)
    df_shrink <- QuickPercentExpression(nuclei_subclass, DGE_results$gene, group1, group2, assay = "SCT", slot = "data")
    DGE_results <- DGE_results |> 
      left_join(df_shrink, by = "gene") |> 
      mutate(log2FoldChange.shrink = log2FoldChange * shrink_coef) |> 
      relocate(gene, log2FoldChange, log2FoldChange.shrink, padj, classification)
    write_csv(DGE_results, glue("{csv_write_dir}/{subclass_fname}__{group1}_vs_{group2}.csv"))
  }
}
