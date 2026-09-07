# Phase 0 of the reviewer-response analysis. Run this first -- every other
# 05/07/08 script reads from the per-subclass
# cache it writes here instead of the full Seurat object.
#
# The full Dec2024 Seurat object is ~17 GB in memory, which leaves no headroom
# on this machine (18 GB RAM) for subset() copies once it's loaded. This script
# loads the object exactly once, strips SCT scale.data (never used downstream),
# and subsets out each subclass to its own small .Rds cache file. The full
# object is then dropped so the later scripts can each run in a small memory
# footprint.

library(Seurat)
library(tidyverse)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/dge_method_comparison"
cache_dir <- SubclassCacheDir(create = TRUE)

nuclei <- LoadDataset("Dec2024")
subclass_list <- LoadSubclassesToUse(nuclei)

# make sure sex is available at the cell level (merge from experimental
# metadata if the object predates the sex column)
if (!"sex" %in% colnames(nuclei@meta.data)) {
  experimental_meta <- read_csv("01-documentation/sample_experimental_metadata_Dec2024.csv",
                                show_col_types = FALSE)
  nuclei$sex <- experimental_meta$sex[match(nuclei$sample, experimental_meta$sample)]
}

# drop the dense scale.data layer; nothing downstream reads it
try(nuclei[["SCT"]]$scale.data <- NULL, silent = TRUE)
gc()

for (subclass in subclass_list) {
  subclass_fname <- ShrinkSubclassName(subclass)
  cache_path <- file.path(cache_dir, glue("{subclass_fname}.Rds"))
  message(glue("Caching {subclass} -> {cache_path}"))

  nuclei_subclass <- subset(nuclei, subclass_name == subclass)
  saveRDS(nuclei_subclass, cache_path)
  rm(nuclei_subclass)
  gc()
}

rm(nuclei)
gc()

print(glue("Phase 0 complete: {length(subclass_list)} subclasses cached in {cache_dir}/"))
