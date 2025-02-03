# load libs & data ----
library(Seurat)
library(tidyverse)
library(reshape2)

source("03-scripts/R/seq_functions.R")

nuclei <- LoadDataset("Dec2024")

csv_write_dir <- "04-analysis/DEGs/Dec2024_activity_condition"


# removeSexSpecificGenes() ----
# removes any sex-specific genes found in the Yao dataset from DEG list
removeSexSpecificGenes <- function(subclass, naive_gene_list) {
  print(glue("Removing sex-specific genes from {subclass}"))
  
  # prep subclass names
  subclass_fname <- gsub(" ", "_", subclass)
  subclass_fname <- gsub("/", "", subclass_fname)
  
  # load in sex-specific genes for this subclass
  fname <- glue("04-analysis/DEGs/Yao2023_sex_DEGs/{subclass_fname}_DEGs.csv")
  sex_gene_list <- read_csv(fname) |> 
    pull(gene)
  
  # remove any genes from naive_gene_list found in sex_gene_list
  scrubbed_gene_list <- naive_gene_list[!naive_gene_list %in% sex_gene_list]
  
  return(scrubbed_gene_list)
}


# findDEG() ----
# uses MAST to find DEGs between conditions within a subclass and removes sex-specific DEGs
findDEG <- function(seurat_obj, subclass, group1, group2, logFC_threshold=0.25) {
  # print report statement
  print(glue("{subclass} between {group1} and {group2}."))

  # skip is number of cells < 30 per group
  if (sum(seurat_obj$activity_condition == group1) < 30 | 
      sum(seurat_obj$activity_condition == group2) < 30) {
    print(glue("Skipping {subclass} between {group1} and {group2} due to low cell count."))

    return(data.frame())
  }
  
  # find DEGs
  DEGs <- FindMarkers(
    seurat_obj, 
    logfc.threshold = logFC_threshold,
    test.use = 'MAST',
    ident.1 = group1, ident.2 = group2,
    latent.vars = c('sublibrary', 'percent.mt'),
  ) |> 
    rownames_to_column(var = 'gene') |> 
    arrange(desc(avg_log2FC))
  
  return(DEGs)
}


# run DEG analysis within each subclasses ----
# select major subclasses ----
subclass_list <- LoadSubclassesToUse(nuclei) |> 
  print(n = 30) |>  # print n before pulling just the list of names
  pull(subclass_name)

log2FC_threshold <- 0.585  # corresponding to 50% change

for (i in seq_along(subclass_list)) {
  subclass <- subclass_list[i]
  print(glue("subclass {i}/{length(subclass_list)}"))

  # subset nuclei before passing to fxn
  nuclei_subclass <- subset(nuclei, subclass_name == subclass)
  Idents(nuclei_subclass) <- nuclei_subclass$activity_condition
  subclass_fname <- ShrinkSubclassName(subclass)

  # get contrasts available for this subclass
  df_contrasts <- GetSubclassContrasts(nuclei, subclass, cell_cutoff = 30)
  
  # run contrasts for cell groups that meet cell_cutoff
  for (i in seq_along(df_contrasts$group1)) {
    group1 = df_contrasts$group1[i]
    group2 = df_contrasts$group2[i]

    DGE_results <- findDEG(nuclei_subclass, subclass, group1, group2, log2FC_threshold)
    write_csv(DGE_results, glue("{csv_write_dir}/{subclass_fname}__{group1}_vs_{group2}.csv"))
  }
}

