# This script plots the distribution of DG Supertypes along anatomical axes.
source("03-scripts/R/seq_functions.R")

library(ggridges)
library(patchwork)
library(ggh4x)
library(gt)

condition_colors <- LoadActivityColors("May2024")
subclass_colors <- LoadAllenColors("subclass")
supertype_colors <- LoadAllenColors("supertype")
cluster_colors <- LoadAllenColors("cluster")
activity_colors <- LoadActivityColors()

supertype_colors_dg <- c(
  "0136 DG Glut_1" = "#009BCC",
  "0137 DG Glut_2" = "#F573FF",
  "0138 DG Glut_3" = "#ff9d2dff",
  "0139 DG Glut_4" = "#c75412ff"
)

cluster_colors_dg <- c(
  "0502 DG Glut_1" = "#972E99",
  "0503 DG Glut_1" = "#5CAECC",
  "0504 DG Glut_1" = "#FF0700",

  "0505 DG Glut_2" = "#CC603D",
  "0506 DG Glut_2" = "#FFD826",
  "0507 DG Glut_2" = "#99176A",

  "0508 DG Glut_3" = "#2E9964",
  "0509 DG Glut_3" = "#CC3DA3",

  "0510 DG Glut_4" = "#8EFF4D"
)


# merge MERFISH meta ----------------------------------------
allen_taxonomy <- read_excel("02-data/published_data/Yao2023/allen_taxonomy_metadata.xlsx")
allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv")
meta_merfish <- read_csv("02-data/published_data/Zhang2023/cell_metadata.csv")

# IMPORTANT:
# "cluster_alias" col in merfish meta is equal to "cl" col in the allen_taxonomy
# "cl" and "cluster_id" are NOT the same thing in the allen_taxonomy. Use "cl" for joins!
allen_taxonomy <- allen_taxonomy |> 
  mutate(cl = as.numeric(cl))
dg_merfish <- meta_merfish |> 
  left_join(allen_taxonomy, by = c("cluster_alias" = "cl")) |> 
  filter(low_quality_mapping == FALSE) |> 
  filter(subclass_id_label == '037 DG Glut') |> 
  filter(x < 9 | x > 2.5) |>      # anything outside this range is mis-classified
  filter(y < 8.1 | y > 2.9 ) |>   # anything outside is mis-classified
  filter(z < 8 & z > 4)           # anything outside is mis-classified

theme_S4 <- theme(
    plot.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none',
    legend.text = element_text(size = 20)
)


# p_DG_supertypes_APaxis ------------------------------------------------
p_DG_supertypes_APaxis <- dg_merfish |> 
  summarise(n = n(), .by = c(z, supertype_id_label)) |> 
  tidyr::complete(z, supertype_id_label, fill = list(n = 0)) |> 
  mutate(composition = n / sum(n), .by = z) |> 
  arrange(z, supertype_id_label) |> 
  mutate(height = composition) |> 
ggplot() +
  geom_area(aes(x = z, y = height, fill = supertype_id_label), position = 'stack') +
  geom_vline(xintercept = seq(7.5, 4.2, -0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = supertype_colors_dg) +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_void() +
  theme_S4
print(p_DG_supertypes_APaxis)


# p_DG_supertypes_DVaxis ------------------------------------------------
# plot DG supertypes along DV axis
p_DG_supertypes_DVaxis <- dg_merfish |> 
  mutate(ycut = cut(y, breaks = seq(2.9, 8.1, 0.2))) |> 
  mutate(yy = as.numeric(substr(as.character(ycut), 2, 4))) |>   # change to number
  summarise(n = n(), .by = c(yy, supertype_id_label)) |>
  tidyr::complete(yy, supertype_id_label, fill = list(n = 0)) |>
  mutate(composition = n / sum(n), .by = yy) |> 
  arrange(yy, supertype_id_label) |> 
  mutate(height = composition) |>
ggplot() +
  geom_area(aes(x = yy, y = height, fill = supertype_id_label), position = 'stack') +
  geom_vline(xintercept = seq(3, 8, 0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = supertype_colors_dg) +
  coord_flip() +
  scale_x_reverse()
  # theme_S4
print(p_DG_supertypes_DVaxis)


# p_DG_clusters_APaxis ------------------------------------------------
p_DG_clusters_APaxis <- dg_merfish |> 
  summarise(n = n(), .by = c(z, cluster_id_label)) |> 
  tidyr::complete(z, cluster_id_label, fill = list(n = 0)) |> 
  mutate(composition = n / sum(n), .by = z) |> 
  arrange(z, cluster_id_label) |> 
  mutate(height = composition) |> 
ggplot() +
  geom_area(aes(x = z, y = height, fill = cluster_id_label), position = 'stack') +
  geom_vline(xintercept = seq(7.5, 4.2, -0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = cluster_colors_dg) +
  scale_x_reverse() +
  scale_y_reverse() +
  theme_void() +
  theme_S4
print(p_DG_clusters_APaxis)


# p_DG_clusters_DVaxis ------------------------------------------------
# plot DG clusters along DV axis
p_DG_clusters_DVaxis <- dg_merfish |> 
  mutate(ycut = cut(y, breaks = seq(2.9, 8.1, 0.2))) |> 
  mutate(yy = as.numeric(substr(as.character(ycut), 2, 4))) |>   # change to number
  summarise(n = n(), .by = c(yy, cluster_id_label)) |>
  tidyr::complete(yy, cluster_id_label, fill = list(n = 0)) |>
  mutate(composition = n / sum(n), .by = yy) |> 
  arrange(yy, cluster_id_label) |> 
  mutate(height = composition) |> 
ggplot() +
  geom_area(aes(x = yy, y = height, fill = cluster_id_label), position = 'stack') +
  geom_vline(xintercept = seq(3, 8, 0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = cluster_colors_dg) +
  coord_flip() +
  scale_x_reverse() +
  theme_void() +
  theme_S4
print(p_DG_clusters_DVaxis)


# if (SAVE_PLOTS == TRUE) {
#   save_plots <- "05-results/LION/raw_R_plots"
#   # png
#   ggsave(plot = p_supertypes, 
#          path = save_plots, 
#          filename = "CA1_supertypes.png",
#          device = png, width = 16, height = 5, dpi = 300)
#   ggsave(plot = p_spatial_AP,
#          path = save_plots,
#          filename = "CA1_supertypes_APaxis.png",
#          device = png, width = 9, height = 1.5, dpi = 300)
#   ggsave(plot = p_spatial_DV,
#          path = save_plots,
#          filename = "CA1_supertypes_DVaxis.png",
#          device = png, width = 16, height = 4, dpi = 300)
#   # svg
#   ggsave(plot = p_supertypes, 
#          path = save_plots, 
#          filename = "CA1_supertypes.svg",
#          width = 9, height = 3, units = 'in')
#   ggsave(plot = p_CA1_supertypes_APaxis,
#          path = save_plots,
#          filename = "CA1_supertypes_APaxis.svg",
#          width = 9, height = 1.5, units = 'in')
#   ggsave(plot = p_CA1_supertypes_DVaxis,
#          path = save_plots,
#          filename = "CA1_supertypes_DVaxis.svg",
#          width = 9, height = 1, units = 'in')
# }
