## ----SuppFig6
# Shows baseline (SE) IEG expression across circadian time (ZT0, ZT4, ZT12, ZT16)
# Addresses r1-comment2: focused visualization of IEG x ZT at baseline
# Per subclass: z-scored heatmap (relative circadian pattern) +
#               bar annotation of absolute normalized counts (as in Figs 6A-D side panels)

library(DESeq2)
library(ComplexHeatmap)
library(circlize)

source('03-scripts/R/seq_functions.R')

if (!exists('nuclei')) { nuclei <- LoadDataset("Dec2024") }
zt_colors <- LoadZTColors()

# Gene and subclass definitions ----------------------------------------
genes_to_plot <- LoadGeneList('IEG')

subclass_list <- c(
  '016 CA1-ProS Glut',
  '037 DG Glut',
  '319 Astro-TE NN',
  '327 Oligo NN'
)

# SE-only subsets
nuclei_ca1 <- nuclei |> subset(subclass_name == subclass_list[1] & activity_condition == 'SE')
nuclei_dg <- nuclei |> subset(subclass_name == subclass_list[2] & activity_condition == 'SE')
nuclei_astro <- nuclei |> subset(subclass_name == subclass_list[3] & activity_condition == 'SE')
nuclei_oligo <- nuclei |> subset(subclass_name == subclass_list[4] & activity_condition == 'SE')

seurat_subsets <- list(
  CA1 = nuclei |> subset(subclass_name == subclass_list[1] & activity_condition == 'SE'),
  DG = nuclei |> subset(subclass_name == subclass_list[2] & activity_condition == 'SE'),
  Astro = nuclei |> subset(subclass_name == subclass_list[3] & activity_condition == 'SE'),
  Oligo = nuclei |> subset(subclass_name == subclass_list[4] & activity_condition == 'SE')
)

# DESeq2 normalization (SE-only, design = ~ ZT) ----------------------------------------
dds_subsets <- list()

for (subclass in names(seurat_subsets)) {
  nuclei_subclass <- seurat_subsets[[subclass]]

  # genes expressed in >1% of cells, plus all IEGs
  mat <- nuclei_subclass[["SCT"]]@data
  expressed_genes <- c(rownames(mat)[rowMeans(mat > 0) > 0.01], genes_to_plot) |> unique()

  pseudobulk_counts <- AggregateExpression(
    nuclei_subclass,
    group.by  = "sample",
    features  = expressed_genes,
    return.seurat = FALSE
  )$RNA

  col_data <- nuclei_subclass@meta.data |>
    as_tibble() |>
    distinct(sample, ZT) |>
    mutate(sample = str_replace(sample, '_', '-')) |>
    column_to_rownames("sample") |>
    arrange(ZT)

  col_data <- col_data[colnames(pseudobulk_counts), , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(
    countData = as.matrix(pseudobulk_counts),
    colData   = col_data,
    design    = ~ ZT
  )
  dds <- DESeq(dds)
  dds_subsets[[subclass]] <- dds
}


# SE IEG heatmaps ----------------------------------------
plots <- list()

for (subclass in names(dds_subsets)) {
  dds      <- dds_subsets[[subclass]]
  col_data <- as.data.frame(colData(dds))

  # normalized counts for IEGs
  norm_counts <- counts(dds, normalized = TRUE) |>
    t() |>
    as.data.frame() |>
    rownames_to_column('sample') |>
    select(sample, all_of(genes_to_plot)) |>
    tibble()

  df_bulk <- col_data |>
    rownames_to_column('sample') |>
    left_join(norm_counts, by = 'sample') |>
    pivot_longer(
      cols      = all_of(genes_to_plot),
      names_to  = 'gene',
      values_to = 'count'
    )

  # mean normalized count per gene × ZT
  gene_means_by_zt <- df_bulk |>
    group_by(gene, ZT) |>
    summarize(mean_count = mean(count, na.rm = TRUE), .groups = 'drop') |>
    pivot_wider(names_from = ZT, values_from = mean_count) |>
    column_to_rownames('gene') |>
    as.matrix()

  # enforce gene order
  gene_means_by_zt <- gene_means_by_zt[genes_to_plot, , drop = FALSE]

  # z-score per gene (row-wise) for heatmap body
  mat_z <- t(scale(t(gene_means_by_zt)))

  # heatmap dimensions
  cell_dim <- 1.5
  hm_h <- unit(nrow(mat_z) * cell_dim, "cm")
  hm_w <- unit(ncol(mat_z) * cell_dim, "cm")

  # ZT column annotation
  zt_levels <- colnames(mat_z)
  col_anno <- HeatmapAnnotation(
    ZT = anno_simple(
      zt_levels,
      col    = setNames(LoadZTColors(), zt_levels),
      height = unit(0.4, "cm")
    ),
    annotation_name_side = 'left',
    annotation_name_gp   = gpar(fontsize = 12)
  )

  # bar annotation — absolute normalized counts (as in Figs 6A-D side panels)
  bar_anno <- rowAnnotation(
    Expression = anno_barplot(
      gene_means_by_zt,
      gp         = gpar(fill = LoadZTColors()),
      width      = unit(cell_dim * 1.5, "cm"),
      axis_param = list(side = "bottom"),
      beside     = TRUE
    )
  )

  p <- Heatmap(
    mat_z,
    name = 'z-score',
    col  = circlize::colorRamp2(c(-2, 0, 2), hcl_palette = 'Blue-Red 3'),
    row_names_gp         = gpar(fontsize = 16, fontface = 'italic'),
    row_names_side       = 'left',
    column_names_gp      = gpar(fontsize = 18),
    column_names_rot     = 45,
    heatmap_legend_param = list(title = 'z-score'),
    cluster_rows         = FALSE,
    cluster_columns      = FALSE,
    column_title         = subclass,
    column_title_gp      = gpar(fontsize = 18),
    bottom_annotation       = col_anno,
    rect_gp              = gpar(col = 'black', lwd = 0.5),
    width                = hm_w,
    height               = hm_h,
    show_heatmap_legend  = TRUE
  ) + bar_anno

  plots <- c(plots, list(p))

  if (SAVE_PLOTS) {
    save_dir <- "05-results/archive/Figure5_SuppFig1_SE_IEG_by_ZT/raw_R_plots"

    save_h <- unit(nrow(mat_z) * cell_dim * 1.15, "cm")
    save_w <- unit(ncol(mat_z) * cell_dim * 2.5,  "cm")

    plot_path <- glue("{save_dir}/SE_IEG__{subclass}.png")
    png(plot_path, width = save_w, height = save_h, units = "cm", res = 900)
    draw(p, heatmap_legend_side = 'right')
    dev.off()

    svgsave(
      plot     = p,
      filename = glue("SE_IEG_{subclass}.svg"),
      path     = save_dir,
      width    = as.numeric(save_w),
      height   = as.numeric(save_h)
    )
  }
}

plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]]
## ----
