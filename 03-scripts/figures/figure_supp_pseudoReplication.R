library(ggplot)
library(patchwork)
library(ggrepel)
library(dplyr)
library(glue)

library(Seurat)


# setup ---------------------------------------------------------------
source('03-scripts/R/seq_functions.R')
subclasses <- LoadSubclassesToUse(nuclei)
nuclei <- LoadDataset('Dec2024')
nuclei <- nuclei |> subset(subclass_name %in% subclasses)


findDEG <- function(seurat_obj, subclass, group1, group2, logFC_threshold=1.2, pct = 0.01) {
  # print report statement
  print(glue("{subclass} between {group1} and {group2}."))
  
  Idents(seurat_obj) <- seurat_obj$activity_condition
  DGE_results <- FindMarkers(
    seurat_obj,
    test.use = 'wilcox',
    ident.1 = group1,
    ident.2 = group2,
    logfc.threshold = logFC_threshold,
    min.pct = pct
  )
  
  DGE_results <- DGE_results |> 
    rownames_to_column('gene') |>
    relocate(gene) |> 
    arrange(desc(avg_log2FC))
  return(DGE_results)
}


# find DEGs with wilcoxson ------------------------------------------------
csv_write_dir <- "04-analysis/DEGs/Dec2024_activity_condition_wilcoxson"
log_FC_threshold <- 0.585

for (subclass in subclasses) {
  subclass_fname <- ShrinkSubclassName(subclass)
  nuclei_subclass <- nuclei |> subset(subclass_name == subclass)
  
  # get gene list to test for this subclass
  mat <- nuclei_subclass[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- rownames(mat)
  
  # get contrasts available for this subclass
  df_contrasts <- GetSubclassContrasts(nuclei, subclass, cell_cutoff = 30)
  
  # run contrasts for cell groups that meet cell_cutoff
  for (k in seq_along(df_contrasts$group1)) {
    group1 <- as.character(df_contrasts$group1[k])
    group2 <- as.character(df_contrasts$group2[k])
    
    print(glue("{subclass} between {group1} and {group2}."))
    
    Idents(nuclei_subclass) <- nuclei_subclass$activity_condition
    DGE_results <- FindMarkers(
      nuclei_subclass,
      features = gene_list,
      test.use = 'wilcox',
      ident.1 = group1,
      ident.2 = group2,
      logfc.threshold = log_FC_threshold,
      min.pct = 1e-20
    )
    
    DGE_results <- DGE_results |> 
      rownames_to_column('gene') |>
      relocate(gene) |> 
      arrange(desc(avg_log2FC))
    
    n_signif <- DGE_results |> 
      filter(p_val_adj < 0.05, abs(avg_log2FC) >= log_FC_threshold) |> 
      nrow()
    print(glue("Number signif genes: {n_signif}"))
    
    write_csv(DGE_results, glue("{csv_write_dir}/{subclass_fname}__{group1}_vs_{group2}.csv"))
  }
}


# load DEGs together ------------------------------------------------
df_n <- tibble(
  subclass = character(),
  n_deg_pseudo = numeric(),
  n_deg_wilcoxson = numeric()
)
for (subclass in subclasses) {
  dir_pseudo <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
  dir_wilcoxson <- "04-analysis/DEGs/Dec2024_activity_condition_wilcoxson"
  subclass_shrink <- ShrinkSubclassName(subclass)
  subclass_fname <- glue('{subclass_shrink}__EE30m_vs_SE.csv')
  
  log_FC_threshold <- 0.585
  n_deg_pseudobulk <- read_csv(glue("{dir_pseudo}/{subclass_fname}"), show_col_types = F) |> 
    filter(padj < 0.05, abs(log2FoldChange) >= log_FC_threshold) |>
    nrow()
  n_deg_wilcoxson  <- read_csv(glue("{dir_wilcoxson}/{subclass_fname}"), show_col_types = F) |> 
    filter(p_val_adj < 0.05, abs(avg_log2FC) >= log_FC_threshold) |>
    nrow()
  
  df_n <- rbind(
    df_n, 
    tibble(subclass = subclass, 
           n_deg_pseudo = n_deg_pseudobulk, 
           n_deg_wilcoxson = n_deg_wilcoxson)
  )
}

# compare number DEGs with each method --------------------------------
df_n |> 
ggplot(aes(n_deg_pseudo, n_deg_wilcoxson)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = 'dashed') +
  labs(x = 'Number DEGs (pseudobulk)', y = 'Number DEGs (wilcoxson)') +
  theme_minimal() +
  theme(legend.position = 'none')