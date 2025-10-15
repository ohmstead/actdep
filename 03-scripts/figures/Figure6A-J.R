## ----Fig6
# Makes hm plots for Clock & IEG expression in SE subclasses
# Makes hm plots for significant interaction genes in CA1 and DG
# Uses DESeq2 to normalize expression values

library(DESeq2)
library(ComplexHeatmap)
library(circlize)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset("Dec2024")
zt_colors <- LoadZTColors()
activity_colors <- LoadActivityColors()

# define genes to test ----------------------------------------
a_priori_genes <- c(LoadGeneList('circadian'), LoadGeneList('IEG'))

subclass_list <- subclass_list <- c(
  '016 CA1-ProS Glut',
  '037 DG Glut',
  '319 Astro-TE NN',
  '327 Oligo NN'
)

nuclei <- LoadDataset('Dec2024')
seurat_subsets <- list(
  CA1   = nuclei |> subset(subclass_name == subclass_list[1]  & activity_condition %in% c('SE', 'EE30m')),
  DG    = nuclei |> subset(subclass_name == subclass_list[2]  & activity_condition %in% c('SE', 'EE30m')),
  Astro = nuclei |> subset(subclass_name == subclass_list[3] & activity_condition %in% c('SE', 'EE30m')),
  Oligo = nuclei |> subset(subclass_name == subclass_list[4] & activity_condition %in% c('SE', 'EE30m'))
)

dds_subsets <- list()

# norm counts with DESeq ---------------------------------------------
for (subclass in names(seurat_subsets)) {
  nuclei_subclass <- seurat_subsets[[subclass]]
  Idents(nuclei_subclass) <- nuclei_subclass$activity_condition
  
  # only test genes expressed in 5% of cells
  mat <- nuclei_subclass[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  expressed_genes <- rownames(mat)
  # ensure all a priori genes are in there too
  expressed_genes <- c(expressed_genes, a_priori_genes) |> unique()
  
  pseudobulk_counts <- AggregateExpression(nuclei_subclass, 
                                           group.by = "sample", 
                                           features = expressed_genes, 
                                           return.seurat = FALSE)$RNA
  col_data <- nuclei_subclass@meta.data |> 
    as_tibble() |> 
    distinct(sample, activity_condition, ZT) |> 
    mutate(sample = str_replace(sample, '_', '-')) |> 
    column_to_rownames("sample") |> 
    arrange(ZT, activity_condition)
  
  col_data$sample <- rownames(col_data)
  col_data <- col_data[colnames(pseudobulk_counts),]
  
  dds <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_counts), 
                                colData = col_data, 
                                design = ~ activity_condition + ZT + activity_condition:ZT)
  
  dds <- DESeq(dds)
  dds_subsets[[subclass]] <- dds
}


# IEG hm ----------------------------------------
genes_to_plot <- LoadGeneList('IEG')
for (subclass in names(dds_subsets)) {
  # get hm count data
  dds <- dds_subsets[[subclass]]
  
  # Get counts
  bulk_counts <- counts(dds, normalized = TRUE) |> 
    t() |> 
    as.data.frame() |>
    rownames_to_column('sample') |> 
    select(sample, all_of(genes_to_plot)) |> 
    tibble()
  
  # 1) build df
  df_bulk <- col_data |> 
    left_join(bulk_counts) |> 
    pivot_longer(
      cols = -c(sample, activity_condition, ZT), 
      names_to = 'gene', 
      values_to = 'count'
    ) |>
    arrange(gene, ZT, activity_condition, sample) |>
    relocate(gene, activity_condition, ZT, sample) |>
    mutate(sample = factor(sample, levels = unique(sample)))
  
  # 2) calculate log2FC
  df_means <- df_bulk |>
    group_by(gene, activity_condition, ZT) |>
    summarize(mean_count = mean(count, na.rm = TRUE) + 1, .groups = "drop") |> 
    group_by(gene, ZT) |> 
    summarize(log2FC = log2(
      mean_count[activity_condition == 'EE30m'] / 
        mean_count[activity_condition == 'SE']
    ))
  
  # 3) prepare matrix for ComplexHeatmap
  mat <- df_means |> 
    pivot_wider(names_from = ZT, values_from = log2FC) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  mat <- mat[genes_to_plot, , drop = FALSE]  # order genes
  
  # 4) ZT annotation
  gene_means_by_zt <- df_bulk |> 
    filter(activity_condition == 'SE') |> 
    group_by(gene, ZT) |> 
    summarize(mean_count = mean(count, na.rm = TRUE), .groups = "drop") |> 
    pivot_wider(names_from = ZT, values_from = mean_count) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  
  # ensure row order matches heatmap
  gene_means_by_zt <- gene_means_by_zt[rownames(mat), , drop = FALSE]
  
  # get dim for hm dimensions
  cell_dim <- 1.5
  hm_h <- unit(nrow(mat) * cell_dim, "cm")
  hm_w <- unit(ncol(mat) * cell_dim, "cm")
  
  # make anno
  bar_anno <- rowAnnotation(
    Expression = anno_barplot(
      gene_means_by_zt,
      gp = gpar(fill = LoadZTColors()),
      width = unit(cell_dim, "cm"),
      axis_param = list(side = "bottom"),
      beside = T,
    )
  )
  
  # 5) plot
  p <- Heatmap(
    mat,
    name = 'log2FC',
    # col = colorRamp2(c(-3, 0, 3), c('blue', 'white', 'red')),
    circlize::colorRamp2(c(-3, 0, 3), hcl_palette = 'RdBu', reverse = T),
    row_names_gp = gpar(fontsize = 16, fontface = 'italic'),
    row_names_side = 'left',
    column_names_gp = gpar(fontsize = 18),
    heatmap_legend_param = list(title = 'log2FC'),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = subclass,
    column_title_gp = gpar(fontsize = 18, fontface = 'bold'),
    column_names_rot = 45,
    rect_gp = gpar(col = 'black', lwd = 0.5),
    width = hm_w,
    height = hm_h,
    show_heatmap_legend = F,
    ) + bar_anno
  print(p)
  
  if (SAVE_PLOTS) {
    save_dir <- "05-results/ORCA/raw_R_plots"
    
    # add 15% to height and width for saving
    save_h <- unit(nrow(mat) * cell_dim * 1.15, "cm")
    save_w <- unit(ncol(mat) * cell_dim * 1.75, "cm")
    
    # png
    plot_path <- glue("{save_dir}/IEG__{subclass}.png")
    png(plot_path, width = save_w, height = save_h, units = "cm", res = 900)
    draw(p, heatmap_legend_side = 'bottom')
    dev.off()
    
    # svg
    svgsave(plot = p,
            filename = glue("IEG_{subclass}.svg"), 
            path = save_dir, 
            width = as.numeric(save_w), height = as.numeric(save_h))
  }
}


# Clock hm ----------------------------------------
genes_to_plot <- c('Per1', 'Per2', 'Cry2', 'Clock', 'Bmal1')
for (subclass in names(dds_subsets)) {
  # get hm count data
  dds <- dds_subsets[[subclass]]
  
  # Get counts
  bulk_counts <- counts(dds, normalized = TRUE) |> 
    t() |> 
    as.data.frame() |>
    rownames_to_column('sample') |> 
    select(sample, all_of(genes_to_plot)) |> 
    tibble()
  
  # 1) build df
  df_bulk <- col_data |> 
    left_join(bulk_counts) |> 
    pivot_longer(
      cols = -c(sample, activity_condition, ZT), 
      names_to = 'gene', 
      values_to = 'count'
    ) |>
    arrange(gene, ZT, activity_condition, sample) |>
    relocate(gene, activity_condition, ZT, sample) |>
    mutate(sample = factor(sample, levels = unique(sample)))
  
  # 2) calculate log2FC
  df_means <- df_bulk |>
    group_by(gene, activity_condition, ZT) |>
    summarize(mean_count = mean(count, na.rm = TRUE) + 1, .groups = "drop") |> 
    group_by(gene, ZT) |> 
    summarize(log2FC = log2(
      mean_count[activity_condition == 'EE30m'] / 
        mean_count[activity_condition == 'SE']
    ))
  
  # 3) prepare matrix for ComplexHeatmap
  mat <- df_means |> 
    pivot_wider(names_from = ZT, values_from = log2FC) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  mat <- mat[genes_to_plot, , drop = FALSE]  # order genes
  
  # 4) ZT annotation
  gene_means_by_zt <- df_bulk |> 
    filter(activity_condition == 'SE') |> 
    group_by(gene, ZT) |> 
    summarize(mean_count = mean(count, na.rm = TRUE), .groups = "drop") |> 
    pivot_wider(names_from = ZT, values_from = mean_count) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  
  # ensure row order matches heatmap
  gene_means_by_zt <- gene_means_by_zt[rownames(mat), , drop = FALSE]
  
  # get dim for hm dimensions
  cell_dim <- 1.5
  hm_h <- unit(nrow(mat) * cell_dim, "cm")
  hm_w <- unit(ncol(mat) * cell_dim, "cm")
  
  # assemble annotation
  bar_anno <- rowAnnotation(
    Expression = anno_barplot(
      gene_means_by_zt,
      gp = gpar(fill = LoadZTColors()),
      width = unit(cell_dim, "cm"),
      axis_param = list(side = "bottom"),
      beside = T,
    )
  )
  
  # 5) plot
  p <- Heatmap(
    mat,
    name = 'log2FC',
    circlize::colorRamp2(c(-2, 0, 2), hcl_palette = 'RdBu', reverse = T),
    row_names_gp = gpar(fontsize = 16, fontface = 'italic'),
    row_names_side = 'left',
    column_names_gp = gpar(fontsize = 18),
    heatmap_legend_param = list(title = 'log2FC'),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = subclass,
    column_title_gp = gpar(fontsize = 18, fontface = 'bold'),
    column_names_rot = 45,
    rect_gp = gpar(col = 'black', lwd = 0.5),
    width = hm_w,
    height = hm_h,
    show_heatmap_legend = F,
  ) + bar_anno
  print(p)
  
  if (SAVE_PLOTS) {
    save_dir <- "05-results/ORCA/raw_R_plots"
    
    # add 15% to height and width for saving
    save_h <- unit(nrow(mat) * cell_dim * 1.15, "cm")
    save_w <- unit(ncol(mat) * cell_dim * 1.75, "cm")
    
    # png
    plot_path <- glue("{save_dir}/Clock__{subclass}.png")
    png(plot_path, width = save_w, height = save_h, units = "cm", res = 900)
    draw(p, heatmap_legend_side = 'bottom')
    dev.off()
    
    # svg
    svgsave(plot = p,
            filename = glue("Clock_{subclass}.svg"), 
            path = save_dir, 
            width = as.numeric(save_w), height = as.numeric(save_h))
  }
}


# Interaction hm ----------------------------------------
gene_lists <- list(
  CA1 = c('Cecr2' ,'Daglb', 'Enoph1', 'Poc1b', 'Zfp420', 'Fcho2', 'Gm5820', 'Itpkb', 'Trmt61a', 'Zfp523'),
  DG = c('Kcnn2', 'Napepld')
)

for (subclass in names(gene_lists)) {
  # get hm count data
  dds <- dds_subsets[[subclass]]
  genes_to_plot <- gene_lists[[subclass]]
  
  # Get counts
  bulk_counts <- counts(dds, normalized = TRUE) |> 
    t() |> 
    as.data.frame() |>
    rownames_to_column('sample') |> 
    select(sample, all_of(genes_to_plot)) |> 
    tibble()
  
  # 1) build df
  df_bulk <- col_data |> 
    left_join(bulk_counts) |> 
    pivot_longer(
      cols = -c(sample, activity_condition, ZT), 
      names_to = 'gene', 
      values_to = 'count'
    ) |>
    arrange(gene, ZT, activity_condition, sample) |>
    relocate(gene, activity_condition, ZT, sample) |>
    mutate(sample = factor(sample, levels = unique(sample)))
  
  # 2) calculate log2FC
  df_means <- df_bulk |>
    group_by(gene, activity_condition, ZT) |>
    summarize(mean_count = mean(count, na.rm = TRUE) + 1, .groups = "drop") |> 
    group_by(gene, ZT) |> 
    summarize(log2FC = log2(
      mean_count[activity_condition == 'EE30m'] / 
        mean_count[activity_condition == 'SE']
    ))
  
  # 3) prepare matrix for ComplexHeatmap
  mat <- df_means |> 
    pivot_wider(names_from = ZT, values_from = log2FC) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  mat <- mat[genes_to_plot, , drop = FALSE]  # order genes
  
  # 4) ZT annotation
  gene_means_by_zt <- df_bulk |> 
    filter(activity_condition == 'SE') |> 
    group_by(gene, ZT) |> 
    summarize(mean_count = mean(count, na.rm = TRUE), .groups = "drop") |> 
    pivot_wider(names_from = ZT, values_from = mean_count) |> 
    column_to_rownames('gene') |> 
    as.matrix()
  
  # ensure row order matches heatmap
  gene_means_by_zt <- gene_means_by_zt[rownames(mat), , drop = FALSE]
  
  # get dim for hm dimensions
  cell_dim <- 1.5
  hm_h <- unit(nrow(mat) * cell_dim, "cm")
  hm_w <- unit(ncol(mat) * cell_dim, "cm")
  
  # assemble annotation
  bar_anno <- rowAnnotation(
    Expression = anno_barplot(
      gene_means_by_zt,
      gp = gpar(fill = LoadZTColors()),
      width = unit(cell_dim, "cm"),
      axis_param = list(side = "bottom"),
      beside = T,
    )
  )
  
  # 5) plot
  p <- Heatmap(
    mat,
    name = 'log2FC',
    circlize::colorRamp2(c(-1.5, 0, 1.5), hcl_palette = 'RdBu', reverse = T),
    row_names_gp = gpar(fontsize = 16, fontface = 'italic'),
    row_names_side = 'left',
    column_names_gp = gpar(fontsize = 18),
    heatmap_legend_param = list(title = 'log2FC'),
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    column_title = subclass,
    column_title_gp = gpar(fontsize = 18, fontface = 'bold'),
    column_names_rot = 45,
    rect_gp = gpar(col = 'black', lwd = 0.5),
    width = hm_w,
    height = hm_h,
    show_heatmap_legend = F,
  ) + bar_anno
  print(p)
  
  if (SAVE_PLOTS) {
    save_dir <- "05-results/ORCA/raw_R_plots"
    
    # add 15% to height and width for saving
    save_h <- unit(nrow(mat) * cell_dim * 1.2, "cm")
    save_w <- unit(ncol(mat) * cell_dim * 1.8, "cm")
    
    # png
    plot_path <- glue("{save_dir}/Significant__{subclass}.png")
    png(plot_path, width = save_w, height = save_h, units = "cm", res = 900)
    draw(p, heatmap_legend_side = 'bottom')
    dev.off()
    
    # svg
    svgsave(plot = p,
            filename = glue("Significant_{subclass}.svg"), 
            path = save_dir, 
            width = as.numeric(save_w), height = as.numeric(save_h))
  }
}
## ----