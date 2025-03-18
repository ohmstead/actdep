# This script plots UMAPs for each gigaclass, colored by activity_condition, subclass, and supertype.
# It also plots, for each subclass, the supertype distribution in each activity_condition.

library(ggplot2)
library(ggvenn)
library(readr)
library(dplyr)
library(glue)
library(patchwork)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset('Dec2024')
activity_colors <- LoadActivityColors()


# gigaclass UMAPs ------------------------------------------------
subclass_colors <- LoadAllenColors('subclass')
supertype_colors <- LoadAllenColors('supertype')
cluster_colors <- LoadAllenColors('cluster')
activity_colors <- LoadActivityColors()

subclass_sets <- LoadSubclassesToUse(nuclei, as_gigaclasses = TRUE)
subclass_sets$excitatory <- subclass_sets$excitatory[-2] # rm CA2
subclass_sets$inhibitory <- subclass_sets$inhibitory[-c(2:6)] # rm low-count interneurons
subclass_sets$glia <- subclass_sets$glia[-1] # rm IMNs

# re-embed nuclei
for (gigatype in names(subclass_sets)) {
  nuclei_gigatype <- nuclei |> 
    subset(subclass_name %in% subclass_sets[[gigatype]]) |> 
    RunUMAP(dims = 1:41)
  p1 <- DimPlot(nuclei_gigatype, group.by = 'subclass_name') +
    scale_color_manual(values = subclass_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p2 <- DimPlot(nuclei_gigatype, group.by = 'activity_condition') +
    scale_color_manual(values = activity_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p3 <- DimPlot(nuclei_gigatype, group.by = 'supertype_name') +
    scale_color_manual(values = supertype_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p4 <- DimPlot(nuclei_gigatype, group.by = 'cluster_name') +
    scale_color_manual(values = cluster_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  
  design <- "
  133
  233
  "
  print(p1+p2+p3 + plot_layout(design = design, widths = c(1,0.75,0.75)))
  si(1500,1000,'png')
}


# make df_cluster_counts ------------------------------------------------
df_cluster_counts <- nuclei@meta.data |> 
  filter(subclass_name %in% LoadSubclassesToUse(nuclei)) |>
  group_by(activity_condition, subclass_name, supertype_name) |>
  count(cluster_name, name = 'cluster_count')


# supertype by condition ------------------------------------------------
# excitatory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$excitatory) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = supertype_colors)
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# inhibitory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$inhibitory) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = supertype_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# glia
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$glia) |>
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = supertype_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())


# cluster by condition ------------------------------------------------
# excitatory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$excitatory) |> 
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# inhibitory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$inhibitory) |> 
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())

# glia
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$glia) |>
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())


# DG clusters by condition ------------------------------------------------
df_cluster_counts |> 
  filter(subclass_name == '037 DG Glut') |>
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  scale_fill_manual(values = cluster_colors) +
  theme_void() +
  theme(legend.position = 'none',
        strip.text = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank())
si(400, 1200)

# tmp color codes for visual cross-ref

# 0136 DG Glut_1 0137 DG Glut_2 0138 DG Glut_3 0139 DG Glut_4 
c("#F573FF",      "#009BCC",      "#FF7726",      "#FFBF26")

# 0503 DG Glut_1 0504 DG Glut_1 0505 DG Glut_2 0506 DG Glut_2 0507 DG Glut_2 0508 DG Glut_3 0509 DG Glut_3 0510 DG Glut_4 
c("#5CAECC",      "#FF0700",      "#CC603D",      "#FFD826",      "#99176A",      "#2E9964",      "#CC3DA3",      "#8EFF4D")

