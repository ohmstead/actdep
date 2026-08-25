## ----Fig4C
# This script is purely for producing screenshots of the CA1 supertypes in MERFISH sections.
# It doesn't save any plots.

source("03-scripts/R/seq_functions.R")

library(dplyr)
library(readr)
library(glue)

library(patchwork)
library(ggplot2)
library(plotly)

condition_colors <- LoadActivityColors("May2024")
subclass_colors <- LoadAllenColors("subclass")
supertype_colors <- LoadAllenColors("supertype")


# load MERFISH and metadata  ----------------------------------------
allen_taxonomy <- read_excel("02-data/published_data/Yao2023/allen_taxonomy_metadata.xlsx")
allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv")
meta_merfish <- read_csv("02-data/published_data/Zhang2023/cell_metadata.csv")


# merge MERFISH and metadata  ----------------------------------------
# IMPORTANT:
# "cluster_alias" col in merfish meta is equal to "cl" col in the allen_taxonomy
# "cl" and "cluster_id" are NOT the same thing in the allen_taxonomy. Use "cl" for joins!
allen_taxonomy <- allen_taxonomy |> 
  mutate(cl = as.numeric(cl))
meta_merfish <- meta_merfish |> 
  left_join(allen_taxonomy, by = c("cluster_alias" = "cl")) |> 
  filter(low_quality_mapping == FALSE)


# make 3-D plot ----------------------------------------
p <- meta_merfish |>
  filter(subclass_id_label == "016 CA1-ProS Glut") |>
  filter(z < 8 & z > 4) |> 
  plot_ly(
    x = ~x, y = ~y, z = ~z,
    color = ~supertype_id_label,
    colors = supertype_colors,
    type = "scatter3d",
    mode = "markers",
    marker = list(size = 5)
  )

# print(p)
## ----