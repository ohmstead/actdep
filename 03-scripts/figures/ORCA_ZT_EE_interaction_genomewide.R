# This script will use DESeq2 to run an interaction analysis on whether genes
# are signficiantly changing in expression at different times of day.

library(ggplot2)
library(dplyr)
library(glue)

library(Seurat)
library(DESeq2)
library(ComplexHeatmap)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset("Dec2024")
zt_colors <- LoadZTColors()
activity_colors <- LoadActivityColors()

# establish subclasses to use ------------------------------------
subclass_list <- LoadSubclassesToUse(nuclei)
subclasses <- subclass_list[c(1, 4, 17, 19)]

seurat_subsets <- list(
  CA1   = nuclei |> subset(subclass_name == subclass_list[1]  & activity_condition %in% c('SE', 'EE30m')),
  DG    = nuclei |> subset(subclass_name == subclass_list[4]  & activity_condition %in% c('SE', 'EE30m')),
  Astro = nuclei |> subset(subclass_name == subclass_list[17] & activity_condition %in% c('SE', 'EE30m')),
  Oligo = nuclei |> subset(subclass_name == subclass_list[19] & activity_condition %in% c('SE', 'EE30m'))
)
fig1_degs <- LoadGeneList('DEGs')

for (subclass in names(seurat_subsets)) {
# make DESeq -----------------------------------------------------
  nuclei_subclass <- seurat_subsets[[subclass]]
  Idents(nuclei_subclass) <- nuclei_subclass$activity_condition
  subclass_fname <- ShrinkSubclassName(subclass)
  
  # only test genes expressed in 5% of cells
  mat <- nuclei_subclass[["SCT"]]@data
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- rownames(mat)
  
  pseudobulk_counts <- AggregateExpression(nuclei_subclass, 
                                           group.by = "sample", 
                                           features = gene_list, 
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
  
  
  # run DESeq --------------------------------------------
  dds <- DESeq(dds)
  FC_thresh <- 0.585
  
  # get DEG
  print(resultsNames(dds))
  
  res1 <- results(dds, name = "activity_conditionEE30m.ZTZT4", tidy = T) |> 
    filter(padj < 0.05) |> 
    mutate(contrast = 'ZT4_vs_ZT0') |> print()
  res2 <- results(dds, name = "activity_conditionEE30m.ZTZT12", tidy = T) |> 
    filter(padj < 0.05) |> 
    mutate(contrast = 'ZT12_vs_ZT0') |> print()
  res3 <- results(dds, name = "activity_conditionEE30m.ZTZT16", tidy = T) |> 
    filter(padj < 0.05) |> 
    mutate(contrast = 'ZT16_vs_ZT0') |> print()
  res4 <- results(dds, 
                  contrast = list(c('activity_conditionEE30m.ZTZT12', 'activity_conditionEE30m.ZTZT4')), 
                  tidy = T) |> 
    filter(padj < 0.05) |> 
    mutate(contrast = 'ZT12_vs_ZT4') |> print()
  res5 <- results(dds, 
                  contrast = list(c('activity_conditionEE30m.ZTZT16', 'activity_conditionEE30m.ZTZT4')), 
                  tidy = T) |> 
    filter(padj < 0.05) |> 
    mutate(contrast = 'ZT16_vs_ZT4') |> print()
  res6 <- results(dds, 
                  contrast = list(c('activity_conditionEE30m.ZTZT16', 'activity_conditionEE30m.ZTZT12')), 
                  tidy = T) |>
    filter(padj < 0.05) |> 
    mutate(contrast = 'ZT16_vs_ZT12') |> print()
  
  # assemble results
  res <- rbind(res1, res2, res3, res4, res5, res6) |> 
    dplyr::rename(gene = row) |> print()
  
  if (nrow(res) == 0) {
    next
  }
  
  for (gene in res$gene) {
    plotCounts(dds, gene, intgroup = c("activity_condition", "ZT"))
  }
  
  
  # test ----------------------------------------
  types <- list(
    CA1 = list(
      type1 = c('Cecr2' ,'Daglb', 'Enoph1', 'Poc1b', 'Zfp420'),
      type2 = c('Fcho2', 'Gm5820', 'Itpkb', 'Trmt61a', 'Zfp523')
    ),
    DG = list(
      type1 = c(),
      type2 = c('Kcnn2', 'Napepld')
    ),
    Astro = list(
      type1 = c(),
      type2 = c()
    ),
    Oligo = list(
      type1 = c(),
      type2 = c()
    )
  )
  types <- types[[subclass]]
  
  bulk_counts <- counts(dds, normalized = T) |> 
    t() |> 
    as.data.frame() |>
    rownames_to_column('sample') |> 
    select(sample, res$gene) |> 
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
    mutate(sample = factor(sample, levels = unique(sample)),
           type = ifelse(gene %in% types$type1, 'type1', 'type2'))
  
  # 2) find mean counts per (gene, activity_condition, ZT)
  df_means <- df |>
    group_by(gene, activity_condition, ZT) |>
    summarize(
      mean_count = mean(count),
      x_min = min(as.numeric(sample)) - 0.4,  # a little left margin
      x_max = max(as.numeric(sample)) + 0.4,  # a little right margin
      sem = sd(count) / sqrt(n()),
      type = dplyr::first(type),
      .groups = "drop"
    )
  
  # 3) plot
  p1 <- ggplot(df) +
    aes(x = sample, y = count, fill = ZT) +
    geom_col(position = 'dodge') +
    scale_fill_manual(values = zt_colors) +
    # lines for each group's mean
    geom_segment(
      data = df_means,
      aes(x = x_min, xend = x_max, y = mean_count, yend = mean_count),
      inherit.aes = FALSE,   # don't use x=sample, y=count from ggplot(df)
      color = "black",
      linewidth = 2
    ) +
    facet_wrap(vars(type, gene), nrow=2, scales = 'free_y') +
    labs(title = glue('Significant ZT:Activity interaction genes: {subclass}'),
         x = 'Pseudobulk samples',
         y = 'Normalized pseudobulk counts') +
    theme(axis.text.x = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5))
  print(p1)
  
  
  # line plot ----------------------------------------
  p2 <- df_means |> 
  ggplot() +
    aes(x = ZT, y = mean_count, color = activity_condition, group = activity_condition) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = mean_count - sem, ymax = mean_count + sem), width = 0.1) +
    facet_wrap(vars(type, gene), nrow=2, scales = "free_y") +
    scale_color_manual(values = activity_colors) +
    labs(title = glue('Significant ZT:Activity interaction genes: {subclass}'),
         x = 'ZT',
         y = 'Normalized pseudobulk counts') +
    theme(axis.text.x = element_text(size = 12, angle = 30, hjust = 1),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5))
  print(p2)
  
  # save plots ----------------------------------------
  if (SAVE_PLOTS) {
  # png
  ggsave(plot = p1,
         path = '05-results/ORCA/raw_R_plots/',
         filename = glue('interacting_genes_genomewide_{subclass_fname}_barplot.png'),
         width =  10, height = 4, dpi = 900, bg = 'white')
  ggsave(plot = p2,
         path = '05-results/ORCA/raw_R_plots/',
         filename = glue('interacting_genes_genomewide_{subclass_fname}_lineplot.png'),
         width =  10, height = 4, dpi = 900, bg = 'white')
  # svg
  ggsave(plot = p1 + LoadBarebonesTheme(ticks = 'y'),
         path = '05-results/ORCA/raw_R_plots/',
         filename = glue('interacting_genes_genomewide_{subclass_fname}_barplot.svg'),
         width =  10, height = 4)
  ggsave(plot = p2 + LoadBarebonesTheme(ticks = 'both'),
         path = '05-results/ORCA/raw_R_plots/',
         filename = glue('interacting_genes_genomewide_{subclass_fname}_lineplot.svg'),
         width =  10, height = 4)
  }
}
