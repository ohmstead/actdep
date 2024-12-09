# ---- load libs & constants ----
print("Loading libraries and constants...")

library(Seurat)
library(dplyr)
library(readr)
library(glue)
library(patchwork)
library(ggplot2)
library(ComplexHeatmap)

source("03-scripts/R/seq_functions.R")

nuclei <- LoadDataset("May2024", "combined")

condition_colors <- LoadConditionColors("May2024")
subclass_colors <- LoadAllSubclassColors(clade = "subclass")
subclass_colors <- LoadAllSubclassColors(clade = "supertype")
subclass_colors <- LoadAllSubclassColors(clade = "cluster")

# ---- load data ----

