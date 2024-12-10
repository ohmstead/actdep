# load libs & data -----
source("03-scripts/R/seq_functions.R")

library(dplyr)
library(readr)
library(glue)

library(patchwork)
library(ggplot2)
library(plotly)

condition_colors <- LoadConditionColors("May2024")
subclass_colors <- LoadAllenColors("subclass")
supertype_colors <- LoadAllenColors("supertype")


# load MERFISH and metadata ----
allen_taxonomy <- read_excel("02-data/published_data/allen_taxonomy_metadata.xlsx")
allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv")
meta_merfish <- read_csv("02-data/published_data/Zhang2023/cell_metadata.csv")


# merge MERFISH and metadata ----
# IMPORTANT:
# "cluster_alias" col in merfish meta is equal to "cl" col in the allen_taxonomy
# "cl" and "cluster_id" are NOT the same thing in the allen_taxonomy. Use "cl" for joins!
allen_taxonomy <- allen_taxonomy |> 
  mutate(cl = as.numeric(cl))
meta_merfish <- meta_merfish |> 
  left_join(allen_taxonomy, by = c("cluster_alias" = "cl")) |> 
  filter(low_quality_mapping == FALSE) |> print()


# plot A-P (z) axis distro ----
p_ap <- meta_merfish |> 
  filter(subclass_id_label == '016 CA1-ProS Glut') |>
  filter(x < 9| x > 2.5) |>   # anything outside this range is mis-classified
  filter(y < 8.1 | y > 2.9 ) |>   # anything outside is mis-classified
  filter(z < 8 & z > 4) |>   # anything outside is mis-classified
  group_by(z, supertype_id_label) |> 
  summarise(n = n()) |>
  mutate(percent = n / sum(n)) |> 
  ggplot() +
  geom_area(aes(x = z, y = percent, fill = supertype_id_label), position = 'fill') +
  geom_vline(xintercept = seq(7.5,4,-0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = supertype_colors) +
  scale_x_reverse() +
  labs(title = "Supertypes along A-P axis",
       x = "<-- anterior / posterior -->", y = "Proportion of cells"
  ) +
  theme(
    plot.title = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.position = c(0.2, 0.25),
  )


# plot D-V (y) axis distro ----
merfish_filtered <- meta_merfish |> 
  filter(subclass_id_label == '016 CA1-ProS Glut') |>
  filter(x < 9| x > 2.5) |>   # anything outside this range is mis-classified
  filter(y < 8.1 | y > 2.9 ) |>   # anything outside is mis-classified
  filter(z < 8 & z > 4) |>   # anything outside is mis-classified
  mutate(ycut = cut(y, breaks = seq(2.9, 8.1, 0.2)))
merfish_filtered <- merfish_filtered |>
  mutate(yy = as.numeric(substr(as.character(ycut), 2, 4)))  # change to number

p_dv <- merfish_filtered |>  
  group_by(yy, supertype_id_label) |> 
  summarise(n = n()) |>
  mutate(percent = n / sum(n)) |> 
  ggplot() +
  geom_area(aes(x = yy, y = percent, fill = supertype_id_label), position = 'fill') +
  geom_vline(xintercept = seq(7.8,3.2,-0.5), color = 'white', linetype = 2, alpha = 0.3) +
  scale_x_reverse() +
  coord_flip() +
  scale_fill_manual(values = supertype_colors) +
  labs(title = "Supertypes along D-V axis",
       x = "<-- ventral / dorsal -->", y = "Proportion of cells"
  ) +
  theme(
    plot.title = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    legend.position = 'none',
  )


p_ap + p_dv + plot_layout(widths = c(2,1))