# load libs & data ----
library(Seurat)
library(tidyverse)
library(DESeq2)

source("03-scripts/R/seq_functions.R")
csv_write_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"

nuclei <- LoadDataset("Dec2024")
log2FC_threshold <- 0.585  # corresponding to 1.5 fold-change (i.e. 50%)
subclass_list <- LoadSubclassesToUse(nuclei)


getConstrastResults <- function(dds, contrast_name, threshold) {
  # get results from Wald test
  res_raw <- results(dds, name = contrast_name)
  
  # filter by shrunken lfc
  res_shrink <- lfcShrink(dds = dds, 
                          res = res_raw, 
                          coef = contrast_name, 
                          type = 'apeglm')|> 
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


# run DEG analysis for each subclass --------------------------------------------
# set up progress bar
pb <- progress_bar$new(
  format = glue("  :subclass_name [:bar] :current/:total (:percent) in :elapsed"),
  total = length(subclass_list), clear = FALSE, width= 60)
pb$tick(0) # init progress bar

for (subclass in subclass_list) {
  # report programm
  pb$tick(tokens = list(subclass_name = subclass))

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
  df_contrasts <- GetSubclassContrasts(nuclei, subclass, cell_cutoff = 30) |> 
    filter(group2 == 'SE')
  
  # run contrasts for cell groups that meet cell_cutoff
  for (k in seq_along(df_contrasts$group1)) {
    group1 <- as.character(df_contrasts$group1[k])
    group2 <- as.character(df_contrasts$group2[k])
    
    contrast_to_test <- glue('activity_condition_{group1}_vs_{group2}')

    DGE_results <- getConstrastResults(dds, contrast = contrast_to_test, log2FC_threshold)
    write_csv(DGE_results, glue("{csv_write_dir}/{subclass_fname}__{group1}_vs_{group2}.csv"))
  }
}
