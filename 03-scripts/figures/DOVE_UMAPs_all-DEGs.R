# This script plots UMAP embeddings computed using the UNION of circadian DEGs
# across all subclasses. Each subclass is embedded in the same DEG feature space.

source('03-scripts/R/seq_functions.R')

# ---- load data ----
nuclei <- LoadDataset('Dec2024')

activity_colors <- LoadActivityColors('Dec2024')
subclass_colors <- LoadAllenColors('subclass')
zt_colors <- LoadZTColors()

# ---- per-subclass circadian DEGs (ZT in SE) ----
deg_dir <- '04-analysis/DEGs/Dec2024_ZT_SE_pseudobulk'

target_subclasses <- c(
  '016 CA1-ProS Glut',
  '037 DG Glut',
  '319 Astro-TE NN',
  '327 Oligo NN'
)

# map subclass names to circadian DEG file stems
group_stem <- c(
  '016 CA1-ProS Glut' = 'CA1',
  '037 DG Glut'       = 'DG',
  '319 Astro-TE NN'   = 'Astro',
  '327 Oligo NN'      = 'Oligo'
)
deg_files <- setNames(file.path(deg_dir, paste0(unname(group_stem[names(group_stem)]), '_all_DEGs.csv')),
                      names(group_stem))

deg_genes_by_subclass <- lapply(names(deg_files), function(s) {
  f <- deg_files[[s]]
  stopifnot(file.exists(f))
  read_csv(f, show_col_types = FALSE) |>
    filter(classification != 'no_change') |>
    pull(gene) |>
    unique()
})
names(deg_genes_by_subclass) <- names(deg_files)

# Build the union of circadian DEGs across all subclasses
union_degs <- sort(unique(unlist(deg_genes_by_subclass)))
union_degs_in_data <- intersect(union_degs, rownames(nuclei))

# ---- per-subclass UMAPs (only that subclass' cells), colored by activity_condition ----
target_subclasses <- c(
  '016 CA1-ProS Glut',
  '037 DG Glut',
  '319 Astro-TE NN',
  '327 Oligo NN'
)

plots <- list()

# Use a fixed number of PCs across subclasses for comparability
n_pcs_global <- min(30, length(union_degs_in_data))

for (sub in target_subclasses) {
  # Create a fresh Seurat object containing only this subclass' cells
  nuclei_sub <- subset(nuclei, subclass_name == sub) |>
    subset(activity_condition == 'SE')
  nuclei_sub <- NormalizeData(nuclei_sub, verbose = FALSE)

  # Recompute scaling and reductions using the UNION of circadian DEGs across subclasses
  common_genes <- intersect(union_degs_in_data, rownames(nuclei_sub))

  nuclei_sub <- ScaleData(nuclei_sub, features = common_genes, verbose = FALSE)
  nuclei_sub <- RunPCA(nuclei_sub, features = common_genes, npcs = n_pcs_global, reduction.name = 'pca.deg', verbose = FALSE)
  nuclei_sub <- RunUMAP(nuclei_sub, reduction = 'pca.deg', dims = 1:n_pcs_global, reduction.name = 'umap.deg', verbose = FALSE, seed.use = 17)

  p <- DimPlot(
    nuclei_sub,
    reduction = 'umap.deg',
    group.by = 'ZT',
    pt.size = 0.1,
    shuffle = TRUE,
    seed = 17,
    raster = FALSE
  ) +
    scale_color_manual(values = zt_colors) +
    theme_minimal() +
    labs(title = glue('{sub} (union circadian DEGs)'), x = '', y = '') +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      legend.position = 'none'
    )

  print(p)
  plots[[sub]] <- p
}

# ---- save ----
if (exists('SAVE_PLOTS') && isTRUE(SAVE_PLOTS)) {
  save_path <- '05-results/NEWT/raw_R_plots'  # reuse existing folder
  for (nm in names(plots)) {
    fname_base <- glue('UMAP_DEG-embedding_condition__{str_replace_all(nm, " ", "_")}')
    ggsave(path = save_path, filename = paste0(fname_base, '.png'), plot = plots[[nm]], width = 9, height = 9, dpi = 900)
    ggsave(path = save_path, filename = paste0(fname_base, '.svg'), plot = plots[[nm]], width = 9, height = 9)
  }
}

print(glue('Script {basename(sys.frame(1)$ofile)} complete!'))
