source('03-scripts/R/seq_functions.R')

library(DESeq2)
library(ComplexHeatmap)

zt_colors <- LoadZTColors()
activity_colors <- LoadActivityColors()

subclass_list <- c(
  '016 CA1-ProS Glut',
  '037 DG Glut',
  '319 Astro-TE NN',
  '327 Oligo NN'
)
seurat_subsets <- list(
  CA1   = nuclei |> subset(subclass_name == subclass_list[1]  & activity_condition %in% c('SE', 'EE30m')),
  DG    = nuclei |> subset(subclass_name == subclass_list[2]  & activity_condition %in% c('SE', 'EE30m')),
  Astro = nuclei |> subset(subclass_name == subclass_list[3] & activity_condition %in% c('SE', 'EE30m')),
  Oligo = nuclei |> subset(subclass_name == subclass_list[4] & activity_condition %in% c('SE', 'EE30m'))
)

for (subclass in names(seurat_subsets)) {
  nuclei_subclass <- seurat_subsets[[subclass]]
  Idents(nuclei_subclass) <- nuclei_subclass$activity_condition
  
  # only test genes expressed in 5% of cells
  mat <- nuclei_subclass[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  expressed_genes <- rownames(mat)
  # ensure a priori genes are in there too
  a_priori_genes <- LoadGeneList('IEG')
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
  
  # line plot
  # genes_to_plot <- c('Per1', 'Per2', 'Clock', 'Bmal1')
  genes_to_plot <- c('Per1')
  
  bulk_counts <- counts(dds, normalized = T) |> 
    t() |> 
    as.data.frame() |>
    rownames_to_column('sample') |> 
    select(sample, all_of(genes_to_plot)) |> 
    tibble() |> 
    print()
  
  # 1) build df
  df <- col_data |> 
    left_join(bulk_counts) |> 
    pivot_longer(
      cols = -c(sample, activity_condition, ZT), 
      names_to = 'gene', 
      values_to = 'count'
    ) |> 
    arrange(gene, activity_condition, ZT, desc(count)) |> 
    mutate(sample = factor(sample, levels = unique(sample)))
  
  # 2) find mean counts per (gene, activity_condition, ZT)
  df_means <- df |>
    group_by(gene, activity_condition, ZT) |>
    summarize(
      mean_count = mean(count),
      x_min = min(as.numeric(sample)) - 0.4,  # a little left margin
      x_max = max(as.numeric(sample)) + 0.4,  # a little right margin
      sem = sd(count) / sqrt(n()),
      .groups = "drop"
    )
  
  # 3) plot
  p2 <- df_means |> 
    ggplot() +
    aes(x = ZT, y = mean_count, color = activity_condition, group = activity_condition) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_count - sem, ymax = mean_count + sem), width = 0.1) +
    facet_wrap(vars(gene), ncol=3, scales = "free_y") +
    scale_color_manual(values = activity_colors) +
    labs(title = glue('{gene_set} genes in {subclass}'),
         y = 'Expression') +
    theme_classic() +
    theme(axis.text.x = element_text(size = 8, angle = 30),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          strip.background = element_blank(),)
  print(p2)
}
