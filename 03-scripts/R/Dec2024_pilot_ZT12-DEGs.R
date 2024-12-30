# need to find DEGs just for the ZT12 timepoint
# load libs and data ----
library(Seurat)
library(SeuratDisk)
library(tidyverse)
library(reshape2)
library(hues)
library(paletteer)
library(patchwork)

source("03-scripts/R/seq_functions.R")

nuclei <- LoadDataset("Dec2024_pilot")

# subset ZT12 ----
nuclei <- nuclei |> 
  subset(
    ZT == 'ZT12' | 
    activity_condition == 'KA30m' | 
    activity_condition == 'KA6h'
  )

# subset only major subclasses ----
subclass_list <- LoadSubclassesToUse(nuclei, cell_cutoff = 120) |> 
  print() |>  # print n before pulling just the list of names
  pull(subclass_name)

nuclei <- nuclei |> 
  subset(subclass_name %in% subclass_list)

# perform DGE ----
csv_write_dir <- "04-analysis/DEGs/Dec2024_pilot_ZT12_activity_condition"
log2FC_threshold <- 0.585  # corresponding to 50% change

for (i in seq_along(subclass_list)) {
  subclass <- subclass_list[i]
  subclass_fname <- ShrinkSubclassName(subclass)
  
  
  print(glue("subclass {i}/{length(subclass_list)}, contrast 1/4"))
  EE30m_vs_SE <- FindDEGs(nuclei, subclass, 'activity_condition', 'EE30m', 'SE', log2FC_threshold)
  write_csv(EE30m_vs_SE, glue("{csv_write_dir}/{subclass_fname}_EE30m_vs_SE.csv"))
  
  print(glue("subclass {i}/{length(subclass_list)}, contrast 2/4"))
  EE6h_vs_SE <- FindDEGs(nuclei, subclass, 'activity_condition', 'EE6h', 'SE', log2FC_threshold)
  write_csv(EE6h_vs_SE, glue("{csv_write_dir}/{subclass_fname}_EE6h_vs_SE.csv"))
  
  print(glue("subclass {i}/{length(subclass_list)}, contrast 3/4"))
  KA30m_vs_SE <- FindDEGs(nuclei, subclass, 'activity_condition', 'KA30m', 'SE', log2FC_threshold)
  write_csv(KA30m_vs_SE, glue("{csv_write_dir}/{subclass_fname}_KA30m_vs_SE.csv"))
  
  print(glue("subclass {i}/{length(subclass_list)}, contrast 4/4"))
  KA6h_vs_SE <- FindDEGs(nuclei, subclass, 'activity_condition', 'KA6h', 'SE', log2FC_threshold)
  write_csv(KA6h_vs_SE, glue("{csv_write_dir}/{subclass_fname}_KA6h_vs_SE.csv"))
}
