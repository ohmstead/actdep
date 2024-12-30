# This script will use the Yao et al. 2023 dataset from the Allen Brain Institute 
# to generate a list of differentially expressed genes between male and female cells.

# load libs & data ----
library(tidyverse)
library(Seurat)

source("03-scripts/R/seq_functions.R")

nuclei <- LoadDataset("Dec2024_pilot")
yao <- LoadDataset("Yao2023")
yao$subclass_name <- yao$subclass


# identify major subclasses ----
major_subclasses <- nuclei@meta.data |> 
  group_by(subclass_name) |>
  summarise(n = n()) |> 
  filter(n >= 150) |> 
  print(n = 100)
  

# perform DGE analysis ----
n <- nrow(major_subclasses)
i <- 1
for (subclass in major_subclasses$subclass_name) {
  print(glue("Performing DGE analysis for {subclass} ({i}/{n})..."))
  
  # need to set sex as Ident
  Idents(yao) <- yao$donor_sex
  
  # filter yao
  yao_subclass <- yao |> subset(subclass_name == subclass)
  
  # find DEGs
  FindMarkers(
    yao_subclass, 
    logfc.threshold = 1,
    test.use = 'wilcox',
    ident.1 = 'M', ident.2 = 'F',
  ) |> 
    filter(p_val_adj < 0.05) |> 
    rownames_to_column(var = 'gene')
  # degs <- FindDEGs(
  #   seurat_obj = yao, 
  #   subclass = subclass, 
  #   ident_var = 'donor_sex', 
  #   group1 = 'M', group2 = 'F', 
  #   logFC_threshold = 1
  # )

  # save to file
  save_dir <- "04-analysis/DEGs/Yao2023_sex_DEGs"
  if (!dir.exists(save_dir)) {dir.create(save_dir)}
  subclass_fname <- ShrinkSubclassName(subclass_name)
  degs_file <- glue("{save_dir}/subclass_fname_sex_DEGs.csv")
  
  i <- i + 1
}
