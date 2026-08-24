## ----SuppFig7
# This script examines the composition of and effects of ZT / activity on cell
# populations important for adult neurogenesis.

# The approach:
#   1. Define related subclasses. Include oligos as a negative-control subclass.
#   2. Re-embed these subclasses in a new UMAP.
#   3. Assign module scores based on marker gene expression. Visualize scores
#      across subclasses (violin grid), ZT, and UMAP.

source('03-scripts/R/seq_functions.R')

if (!exists('nuclei')) { nuclei <- LoadDataset("Dec2024") }

subclass_colors <- LoadAllenColors('subclass')
supertype_colors <- LoadAllenColors('supertype')
cluster_colors <- LoadAllenColors('cluster')
zt_colors       <- LoadZTColors()

SAVE_PLOTS <- FALSE
save_dir   <- "05-results/SuppFigure7/raw_R_plots"


# 1. Define related subclasses. Include oligos as a negative-control subclass ----
subclasses_to_use <- c("037 DG Glut", "038 DG-PIR Ex IMN", "319 Astro-TE NN", "016 CA1-ProS Glut")

# Neurogenic subclasses ordered stem → mature; Oligo = negative control (UMAP only)
neurogenic_subclasses <- c("319 Astro-TE NN", "038 DG-PIR Ex IMN", "037 DG Glut", "016 CA1-ProS Glut")


# 2. Subset and re-embed in new UMAP ----
nuclei_gc <- nuclei |>
  subset((subclass_name %in% subclasses_to_use) & activity_condition == 'SE') |>
  RunUMAP(dims = 1:30, seed.use = 17)

p_umap_subclass <- DimPlot(
  nuclei_gc, group.by = 'subclass_name',
  shuffle = TRUE, seed = 17, pt.size = 0.3
) +
  scale_color_manual(values = subclass_colors) +
  theme_void() +
  theme(plot.title = element_blank(), legend.position = 'right')


# 3a. Assign module scores ----
# rgl_genes         <- c("Hopx", "Sox2", "Fabp7", "Slc1a3", "Nes", "Aqp4")
rgl_genes         <- c("Hopx", "Prom1", "Notch2", "Hes1", "Lfng", "Rbpj")
ipc_genes         <- c("Ascl1", "Eomes", "Mki67", "Top2a", "Pcna", "Ccnd2")
immature_gc_genes <- c("Dcx", "Neurod1", "Calb2", "Sox4", "Sox11", "Stmn1", "Tubb3", "Igfbpl1")
mature_gc_genes   <- c("Prox1", "Calb1")

nuclei_gc <- nuclei_gc |>
  AddModuleScore(features = list(rgl_genes),         name = "Score_RGL", nbin = 12) |>
  AddModuleScore(features = list(ipc_genes),         name = "Score_IPC", nbin = 12) |>
  AddModuleScore(features = list(immature_gc_genes), name = "Score_iGC", nbin = 12) |>
  AddModuleScore(features = list(mature_gc_genes),   name = "Score_mGC", nbin = 12)

# AddModuleScore appends "1" to the name argument
module_cols   <- c("Score_RGL1", "Score_IPC1", "Score_iGC1", "Score_mGC1")
module_labels <- c(Score_RGL1 = "RGL score", Score_IPC1 = "IPC score",
                   Score_iGC1 = "iGC score", Score_mGC1 = "mGC score")


# 3b. Visualize scores across subclasses, ZT, and UMAP ----

# --- Violin grid: 4 rows (modules) x 3 cols (neurogenic subclasses) ---
# facet_grid with scales = 'free_y' shares y axes within each row (module)
df_meta <- nuclei_gc@meta.data |>
  rownames_to_column('cell') |>
  select(cell, ends_with('_name'), all_of(module_cols)) |>
  pivot_longer(cols = all_of(module_cols), names_to = 'module', values_to = 'score') |>
  mutate(
    module        = factor(module, levels = module_cols, labels = unname(module_labels)),
    subclass_name = factor(subclass_name, levels = neurogenic_subclasses),
    supertype_name = factor(supertype_name),
    cluster_name = factor(cluster_name)
  )

p_violin_grid <- ggplot(df_meta, aes(x = factor(1), y = score, fill = subclass_name)) +
  geom_violin(color = NA, alpha = 0.85) +
  # geom_jitter(width = 0.25, size = 0.2, alpha = 0.15, color = 'gray20') +
  geom_hline(yintercept = 0, linetype = 'dashed', color = 'gray60', linewidth = 0.4) +
  # scale_fill_manual(values = supertype_colors) +
  scale_fill_manual(values = subclass_colors) +
  scale_x_discrete(expand = expansion(add = 0)) +
  facet_grid(rows = vars(module), cols = vars(subclass_name), scales = 'free_y') +
  labs(x = NULL, y = 'module score', size = 12) +
  theme_minimal(base_family = 'Helvetica') +
  theme(
    legend.position    = 'none',
    axis.text.x        = element_blank(),
    axis.ticks.x       = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.spacing.x    = unit(10, 'pt'),
    strip.text.y  = element_text(angle = 270, hjust = 0.5, size = 14),
    strip.text.x  = element_text(angle = 0, hjust = 0.5, size = 14)
  )
p_violin_grid


# --- UMAP colored by module score ---
umap_plots <- lapply(module_cols, function(mod) {
  FeaturePlot(nuclei_gc, features = mod, order = TRUE, pt.size = 0.3) +
    scale_color_gradient(low = 'gray80', high = 'red') +
    labs(title = module_labels[[mod]]) +
    theme_void() +
    theme(plot.title = element_text(size = 11, hjust = 0.5), legend.position = 'none')
})
names(umap_plots) <- module_cols

p_umap_all <- p_umap_subclass |
  umap_plots[[1]] | umap_plots[[2]] | umap_plots[[3]] | umap_plots[[4]]


# Print ----
p_umap_all
p_violin_grid


# Save ----
if (SAVE_PLOTS) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  ggsave(glue("{save_dir}/umap_subclass_and_modules.png"),
         plot = p_umap_all, width = 20, height = 4, dpi = 300, bg = 'white')
  ggsave(glue("{save_dir}/umap_subclass_and_modules.svg"),
         plot = p_umap_all, width = 20, height = 4, bg = 'transparent')

  ggsave(glue("{save_dir}/violin_module_x_subclass.png"),
         plot = p_violin_grid, width = 9, height = 12, dpi = 300, bg = 'white')
  ggsave(glue("{save_dir}/violin_module_x_subclass.svg"),
         plot = p_violin_grid, width = 9, height = 12, bg = 'transparent')
}
## ----
