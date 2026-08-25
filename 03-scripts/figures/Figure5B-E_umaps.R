## ----Fig5B-E_umap
# This script embeds ALL selected nuclei together in a single Seurat object
# using the UNION of circadian DEGs across subclasses (ZT in SE pseudobulk).
# Outputs six UMAPs: one colored by subclass, one by ZT.
# Other 4 UMAPs highlight each subclass indiv (colored by ZT; others grey).

source('03-scripts/R/seq_functions.R')

# ---- load data ----
if (!exists('nuclei')) {nuclei <- LoadDataset('Dec2024')}

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

# ---- single-object UMAP using union circadian DEGs ----
# Combine all target subclasses and restrict to SE (ZT is only defined in SE)
nuclei_combined <- nuclei |>
  subset(subclass_name %in% target_subclasses) |>
  subset(activity_condition == 'SE')

# Only keep DEGs present in this combined object
union_degs_in_data <- intersect(union_degs, rownames(nuclei_combined))

# Normalize, scale, and run reductions on shared circadian DEG space
nuclei_combined <- NormalizeData(nuclei_combined, verbose = FALSE)
nuclei_combined <- ScaleData(nuclei_combined, features = union_degs_in_data, verbose = FALSE)

n_pcs <- min(30, length(union_degs_in_data))
nuclei_combined <- RunPCA(
  nuclei_combined,
  features = union_degs_in_data,
  npcs = n_pcs,
  reduction.name = 'pca.deg',
  seed.use = 17
)
nuclei_combined <- RunUMAP(
  nuclei_combined,
  reduction = 'pca.deg',
  dims = 1:n_pcs,
  reduction.name = 'umap.deg',
  seed.use = 17
)

# Prepare color vectors for present levels only
present_subclasses <- sort(unique(nuclei_combined$subclass_name))
present_zt <- sort(unique(nuclei_combined$ZT))
subclass_colors_present <- subclass_colors[present_subclasses]
zt_colors_present <- zt_colors[present_zt]

# Build plots
plots_omnibus <- list()

p_by_subclass <- DimPlot(
  nuclei_combined,
  reduction = 'umap.deg',
  group.by = 'subclass_name',
  pt.size = 0.1,
  shuffle = TRUE,
  seed = 17,
  raster = FALSE
) +
  scale_color_manual(values = subclass_colors_present, drop = FALSE) +
  theme_minimal() +
  labs(title = 'Subclasses', x = '', y = '') +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    # legend.position = 'none'
  )

p_by_zt <- DimPlot(
  nuclei_combined,
  reduction = 'umap.deg',
  group.by = 'ZT',
  pt.size = 0.1,
  shuffle = TRUE,
  seed = 17,
  raster = FALSE
) +
  scale_color_manual(values = zt_colors_present, drop = FALSE) +
  theme_minimal() +
  labs(title = 'ZT', x = '', y = '') +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    # legend.position = 'none'
  )

# print(p_by_subclass)
# print(p_by_zt)
plots_omnibus[['by_subclass']] <- p_by_subclass
plots_omnibus[['by_ZT']] <- p_by_zt


# ---- plot UMAPs
plots <- list()
for (sub in target_subclasses) {
  # create a masked metadata column where non-target cells are labeled 'Other'.
  # this is to non-target subclass grey in the plot
  nuclei_combined$ZT_masked <- ifelse(
    nuclei_combined$subclass_name == sub,
    as.character(nuclei_combined$ZT),
    'Other'
  )
  nuclei_combined$ZT_masked <- factor(nuclei_combined$ZT_masked)

  # Build color vector: grey for 'Other' + ZT palette for present levels
  present_lvls <- levels(nuclei_combined$ZT_masked)
  zt_lvls <- intersect(names(zt_colors), present_lvls)
  color_map <- c('Other' = 'grey80')
  if (length(zt_lvls) > 0) {
    color_map <- c(color_map, zt_colors[zt_lvls])
  }

  p_masked <- DimPlot(
    nuclei_combined,
    reduction = 'umap.deg',
    group.by = 'ZT_masked',
    pt.size = 0.1,
    shuffle = TRUE,
    seed = 17,
    raster = FALSE
  ) +
    scale_color_manual(values = color_map, drop = FALSE) +
    theme_minimal() +
    labs(title = group_stem[[sub]], x = '', y = '') +
    theme(
      panel.grid = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(hjust = 0.5)
    )

  plots[[glue('highlight__{sub}')]] <- p_masked
}

p <- Reduce(`+`, plots) + plot_layout(ncol = 2, guides = 'collect')
p

# ---- save ----
if (exists('SAVE_PLOTS') & SAVE_PLOTS==T) {
  save_path <- '05-results/Figure5/raw_R_plots'  # reuse existing folder
  ggsave(path = save_path, filename = 'UMAP_DEG-embedding_omnibus__by-subclass.png', plot = plots_omnibus[['by_subclass']], width = 9, height = 9, dpi = 900)
  ggsave(path = save_path, filename = 'UMAP_DEG-embedding_omnibus__by-subclass.svg', plot = plots_omnibus[['by_subclass']], width = 9, height = 9)
  ggsave(path = save_path, filename = 'UMAP_DEG-embedding_omnibus__by-ZT.png', plot = plots_omnibus[['by_ZT']], width = 9, height = 9, dpi = 900)
  ggsave(path = save_path, filename = 'UMAP_DEG-embedding_omnibus__by-ZT.svg', plot = plots_omnibus[['by_ZT']], width = 9, height = 9)

  # Save subclass-highlighted versions
  for (sub in target_subclasses) {
    key <- glue('highlight__{sub}')
    fname_base <- glue('UMAP_DEG-embedding_omnibus__subclass-highlight__{str_replace_all(sub, " ", "_")}')
    ggsave(path = save_path, filename = paste0(fname_base, '.png'), plot = plots[[key]], width = 9, height = 9, dpi = 900)
    ggsave(path = save_path, filename = paste0(fname_base, '.svg'), plot = plots[[key]], width = 9, height = 9)
  }
}
## ----