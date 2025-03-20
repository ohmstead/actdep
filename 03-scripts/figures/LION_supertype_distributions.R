# load libs & data -----
source("03-scripts/R/seq_functions.R")

library(dplyr)
library(readr)
library(readxl)
library(glue)

library(patchwork)
library(ggplot2)
library(plotly)

source("03-scripts/R/seq_functions.R")

condition_colors <- LoadActivityColors("May2024")
subclass_colors <- LoadAllenColors("subclass")
supertype_colors <- LoadAllenColors("supertype")


# load MERFISH and metadata ----
allen_taxonomy <- read_excel("02-data/published_data/Yao2023/allen_taxonomy_metadata.xlsx")
allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv")
meta_merfish <- read_csv("02-data/published_data/Zhang2023/cell_metadata.csv")


# merge MERFISH and metadata ----
# IMPORTANT:
# "cluster_alias" col in merfish meta is equal to "cl" col in the allen_taxonomy
# "cl" and "cluster_id" are NOT the same thing in the allen_taxonomy. Use "cl" for joins!
allen_taxonomy <- allen_taxonomy |> 
  mutate(cl = as.numeric(cl))
ca1_merfish <- meta_merfish |> 
  left_join(allen_taxonomy, by = c("cluster_alias" = "cl")) |> 
  filter(low_quality_mapping == FALSE) |> 
  filter(subclass_id_label == '016 CA1-ProS Glut') |>
  filter(x < 9 | x > 2.5) |>      # anything outside this range is mis-classified
  filter(y < 8.1 | y > 2.9 ) |>   # anything outside is mis-classified
  filter(z < 8 & z > 4)           # anything outside is mis-classified


# plot A-P (z) axis distro ----
p_ap <- ca1_merfish |> 
  group_by(z, supertype_id_label) |> 
  summarise(n = n()) |>
  mutate(percent = n / sum(n)) |> 
ggplot() +
  geom_area(aes(x = z, y = percent, fill = supertype_id_label), position = 'fill') +
  geom_vline(xintercept = seq(7.5,4.2,-0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = supertype_colors) +
  scale_x_reverse() +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none',
  )


# plot D-V (y) axis distro ----
p_dv <- ca1_merfish |>
  mutate(ycut = cut(y, breaks = seq(2.9, 8.1, 0.2))) |> 
  mutate(yy = as.numeric(substr(as.character(ycut), 2, 4))) |>   # change to number
  group_by(yy, supertype_id_label) |> 
  summarise(n = n()) |>
  mutate(percent = n / sum(n)) |> 
ggplot() +
  geom_area(aes(x = yy, y = percent, fill = supertype_id_label), position = 'fill') +
  geom_vline(xintercept = seq(7.8,3.2,-0.5), color = 'white', linetype = 2, alpha = 0.3) +
  scale_x_reverse() +
  coord_flip() +
  scale_fill_manual(values = supertype_colors) +
  theme_void() +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none',
  )

print(p_ap + p_dv + plot_layout(widths = c(2,1)))

if (SAVE_PLOTS) {
  save_path <- "05-results/LION/raw_R_plots"
  # png
  ggsave(plot = p_ap, 
         path = save_path, 
         filename = "AP_supertype_distibution.png",
         device = png, width = 6, height = 3.2, dpi = 900)
  ggsave(plot = p_dv, 
         path = save_path, 
         filename = "DV_supertype_distibution.png",
         device = png, width = 4, height = 3.2, dpi = 900)
  # svg
  ggsave(plot = p_ap, 
         path = save_path, 
         filename = "AP_supertype_distibution.svg",
         width = 6, height = 3.2)
  ggsave(plot = p_dv, 
         path = save_path, 
         filename = "DV_supertype_distibution.svg",
         width = 4, height = 3.2)
}
