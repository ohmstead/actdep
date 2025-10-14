## ----Fig 2C-D
# This script plots UMAPs for each gigaclass, colored by activity_condition, subclass, and supertype.
# It also plots, for each subclass, the supertype distribution in each activity_condition.
source('03-scripts/R/seq_functions.R')

library(ggbreak)

save_dir <- '05-results/Figure2/raw_R_plots/'
if (!exists('nuclei')) {
  nuclei <- LoadDataset('Dec2024')
}
seurat_subsets <- LoadDataset('Dec2024', as_gigaclasses = T)

activity_colors <- LoadActivityColors()
subclass_colors <- LoadAllenColors('subclass')
supertype_colors <- LoadAllenColors('supertype')
cluster_colors <- LoadAllenColors('cluster')

# gigaclass UMAPs ------------------------------------------------
# remove low-N subclasses
subclass_sets <- LoadSubclassesToUse(nuclei, as_gigaclasses = TRUE)
subclass_sets$excitatory <- subclass_sets$excitatory[-2] # rm CA2

# re-embed nuclei
nuclei_gigaclass <- seurat_subsets[['excitatory']] |> 
  RunUMAP(dims = 1:41, seed.use = 17)
p1 <- DimPlot(nuclei_gigaclass, group.by = 'subclass_name', shuffle = T, seed = 17) +
  scale_color_manual(values = subclass_colors) +
  theme_void() +
  theme(legend.position = 'none',
        # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_blank())
p2 <- DimPlot(nuclei_gigaclass, group.by = 'activity_condition', shuffle = T, seed = 17) +
  scale_color_manual(values = activity_colors) +
  theme_void() +
  theme(legend.position = 'none',
        # panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        plot.title = element_blank())
p3 <- DimPlot(nuclei_gigaclass, group.by = 'supertype_name', shuffle = T, seed = 17) +
  scale_color_manual(values = supertype_colors) +
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


# make df_cluster_counts ------------------------------------------------
df_cluster_counts <- nuclei@meta.data |> 
  filter(subclass_name %in% LoadSubclassesToUse(nuclei)) |> 
  group_by(activity_condition, subclass_name, supertype_name) |> 
  dplyr::count(cluster_name, name = 'cluster_count') |> 
  group_by(activity_condition, subclass_name)


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
    axis.text.x = element_text(angle=30,hjust=1),
    legend.position = 'none'
  )
print(p1)

# cluster by condition ------------------------------------------------
# excitatory
p2 <- df_cluster_counts |> 
  filter(subclass_name == '037 DG Glut') |> 
ggplot() +
  aes(x = activity_condition, y = cluster_count, fill = cluster_name) +
  geom_bar(stat = "identity", position = 'fill', linewidth = 1) +
  facet_wrap(~subclass_name, nrow = 1) +
  scale_fill_manual(values = cluster_colors) +
  labs(x = '', y = 'Cluster proportion') +
  theme(legend.position = 'none')
print(p2)


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
## ----