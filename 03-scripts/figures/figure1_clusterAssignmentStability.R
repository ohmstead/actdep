library(tidyverse)
library(reticulate)
library(scCustomize)

reticulate::use_condaenv("sc")
source("03-scripts/R/seq_functions.R")

# load in seurat
nuclei_with <- LoadDataset("Dec2024")

# 1. Run Dec2024 with ARGs through MapMyCells.
      # already done

# 2. Remove Tyssowski ARGs from Seurat object.
genes_to_remove <- LoadGeneList("tyssowski")
features_with <- Features(nuclei_with)
features_without <- features_with[!features_with %in% genes_to_remove]
nuclei_without <- subset(nuclei_with, features = features_without)

# 3. Run Dec2024 without ARGs through MapMyCells.
working_save_dir <- "04-analysis/cluster_reassignment"
as.anndata(x = nuclei_without, file_path = working_save_dir, file_name = "Dec2024_pilot_withoutARGs", main_layer = "counts", other_layers = NULL)

# 4. Merge metadata from both results into a unified dataframe.
# 5. Find all discordant labels at the class, subclass, supertype, and cluster level.
# 6. For each concordant barcode, calculate the difference in bootstrapped cluster-assignment probability with and without IEGs.
