# This script plots UMAPs for each gigaclass, colored by activity_condition, subclass, and supertype.
# It also plots, for each subclass, the supertype distribution in each activity_condition.
source('03-scripts/R/seq_functions.R')

library(ggbreak)

save_dir <- '05-results/CRAB/raw_R_plots/'
nuclei <- LoadDataset('Dec2024')
seurat_subsets <- LoadDataset('Dec2024', as_gigaclasses = T)

activity_colors <- LoadActivityColors()
subclass_colors <- LoadAllenColors('subclass')
supertype_colors <- LoadAllenColors('supertype')
cluster_colors <- LoadAllenColors('cluster')

# gigaclass UMAPs ------------------------------------------------
# remove low-N subclasses
subclass_sets <- LoadSubclassesToUse(nuclei, as_gigaclasses = TRUE)
subclass_sets$excitatory <- subclass_sets$excitatory[-2] # rm CA2
subclass_sets$inhibitory <- subclass_sets$inhibitory[-c(2:6)] # rm low-count interneurons
subclass_sets$glia <- subclass_sets$glia[-1] # rm IMNs

# re-embed nuclei
for (gigaclass in names(seurat_subsets)) {
  nuclei_gigaclass <- seurat_subsets[[gigaclass]] |> 
    RunUMAP(dims = 1:41, seed.use = 17)
  p1 <- DimPlot(nuclei_gigaclass, group.by = 'subclass_name') +
    scale_color_manual(values = subclass_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p2 <- DimPlot(nuclei_gigaclass, group.by = 'activity_condition') +
    scale_color_manual(values = activity_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p3 <- DimPlot(nuclei_gigaclass, group.by = 'supertype_name') +
    scale_color_manual(values = supertype_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  p4 <- DimPlot(nuclei_gigaclass, group.by = 'cluster_name') +
    scale_color_manual(values = cluster_colors) +
    theme_void() +
    theme(legend.position = 'none',
          # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
          plot.title = element_blank())
  
  design <- "
  133
  233
  "
  p_UMAPs <- p1+p2+p3 + plot_layout(design = design, widths = c(1,0.75,0.75))
  print(p_UMAPs)
  
  if (SAVE_PLOTS) {
    # png
    ggsave(plot = p_UMAPs,
           filename = glue('reembedded_UMAPs_{gigaclass}.png'),
           path = save_dir,
           width = 15, height = 10, dpi = 300, bg = 'white')
    # svg
    ggsave(plot = p_UMAPs,
           filename = glue('reembedded_UMAPs_{gigaclass}.svg'),
           path = save_dir,
           width = 15, height = 10, bg = 'transparent')
    
  }
}


# make df_cluster_counts ------------------------------------------------
df_cluster_counts <- nuclei@meta.data |> 
  filter(subclass_name %in% LoadSubclassesToUse(nuclei)) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  dplyr::count(cluster_name, name = 'cluster_count') |> 
  group_by(activity_condition, subclass_name) |> 
  print()


# supertype by condition ------------------------------------------------
# excitatory
p1 <- df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$excitatory) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  scale_fill_manual(values = supertype_colors) +
  labs(x = '', y = 'Supertype proportion') +
  facet_wrap(~subclass_name, nrow = 1) +
  theme(
    axis.text.x = element_text(angle=30,hjust=1)
  )
print(p1)
ggplotly(p1)

# inhibitory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$inhibitory) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  summarize(supertype_count = sum(cluster_count)) |> 
ggplot() +
  aes(x = activity_condition, y = supertype_count, fill = supertype_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = supertype_colors)

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
p2 <- df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$excitatory) |> 
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors)
print(p2)
ggplotly(p2)

# inhibitory
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$inhibitory) |> 
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors)

# glia
df_cluster_counts |> 
  filter(subclass_name %in% subclass_sets$glia) |>
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors)


# DG specifically ------------------------------------------------
p3 <- df_cluster_counts |> 
  ungroup() |> 
  filter(subclass_name == '037 DG Glut') |> 
  mutate(cluster_pct = cluster_count / sum(cluster_count) * 100, 
         .by = c(activity_condition, subclass_name)) |>
  filter(supertype_name == '0138 DG Glut_3') |> 
ggplot() +
  aes(x = activity_condition, y = cluster_pct, fill = cluster_name) +
  geom_col(position = 'stack', linewidth = 1) +
  scale_y_break(c(2, 88)) +
  scale_fill_manual(values = cluster_colors)
print(p3)


# save ----------------------------------------
if (SAVE_PLOTS) {
  # png
  ggsave(plot = p1,
         filename = 'excitatory_supertype_by_condition.png',
         path = save_dir,
         width = 15, height = 10, dpi = 300, bg = 'white')
  ggsave(plot = p2,
         filename = 'excitatory_cluster_by_condition.png',
         path = save_dir,
         width = 15, height = 10, dpi = 300, bg = 'white')
  ggsave(plot = p3,
         filename = 'DG_Glut3_cluster_by_condition.png',
         path = save_dir,
         width = 15, height = 10, dpi = 300, bg = 'white')
  # svg
  ggsave(plot = p1 + LoadBarebonesTheme(ticks = 'y'),
         filename = 'excitatory_supertype_by_condition.svg',
         path = save_dir,
         width = 15, height = 10, bg = 'transparent')
  ggsave(plot = p2 + LoadBarebonesTheme(ticks = 'y'),
         filename = 'excitatory_cluster_by_condition.svg',
         path = save_dir,
         width = 15, height = 10, bg = 'transparent')
  ggsave(plot = p3 + LoadBarebonesTheme(ticks = 'y'),
         filename = 'DG_Glut3_cluster_by_condition.svg',
         path = save_dir,
         width = 15, height = 10, bg = 'transparent')
  
}


# tmp color codes for visual cross-ref
# 0136 DG Glut_1 0137 DG Glut_2 0138 DG Glut_3 0139 DG Glut_4 
c("#F573FF",      "#009BCC",      "#FF7726",      "#FFBF26")

# 0503 DG Glut_1 0504 DG Glut_1 0505 DG Glut_2 0506 DG Glut_2 0507 DG Glut_2 0508 DG Glut_3 0509 DG Glut_3 0510 DG Glut_4 
c("#5CAECC",      "#FF0700",      "#CC603D",      "#FFD826",      "#99176A",      "#2E9964",      "#CC3DA3",      "#8EFF4D")

