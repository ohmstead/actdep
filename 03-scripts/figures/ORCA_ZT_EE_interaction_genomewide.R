# This script will use DESeq2 to run an interaction analysis on whether genes
# are signficiantly changing in expression at different times of day.

source('03-scripts/R/seq_functions.R')

library(DESeq2)
library(ComplexHeatmap)

zt_colors <- LoadZTColors()
activity_colors <- LoadActivityColors()

# establish subclasses to use ------------------------------------
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
  
  res1 <- results(dds, name = "activity_conditionEE30m.ZTZT4") %>%
    as.data.frame() |> 
    as_tibble(rownames = 'gene') |> 
    mutate(contrast = 'ZT4_vs_ZT0')
  res2 <- results(dds, name = "activity_conditionEE30m.ZTZT12") %>%
    as.data.frame() |> 
    as_tibble(rownames = 'gene') |> 
    mutate(contrast = 'ZT12_vs_ZT0')
  res3 <- results(dds, name = "activity_conditionEE30m.ZTZT16") %>%
    as.data.frame() |> 
    as_tibble(rownames = 'gene') |> 
    mutate(contrast = 'ZT16_vs_ZT0')
  res4 <- results(
    dds,
    contrast = list(c('activity_conditionEE30m.ZTZT12', 'activity_conditionEE30m.ZTZT4'))
    ) |> 
    as.data.frame() |> 
    as_tibble(rownames = 'gene') |> 
    mutate(contrast = 'ZT12_vs_ZT4')
  res5 <- results(
    dds,
    contrast = list(c('activity_conditionEE30m.ZTZT16', 'activity_conditionEE30m.ZTZT4'))
    ) |>
    as.data.frame() |> 
    as_tibble(rownames = 'gene') |> 
    mutate(contrast = 'ZT16_vs_ZT4')
  res6 <- results(
    dds,
    contrast = list(c('activity_conditionEE30m.ZTZT16', 'activity_conditionEE30m.ZTZT12'))
    ) |>
    as.data.frame() |> 
    as_tibble(rownames = 'gene') |> 
    mutate(contrast = 'ZT16_vs_ZT12')
  
  # assemble results
  res <- rbind(res1, res2, res3, res4, res5, res6) |> 
    arrange(padj) |> 
    print()
  
  # save as csv
  csv_dir <- "04-analysis/DEGs/Dec2024_ZT-EE_interaction"
  write_csv(res, glue("{csv_dir}/a_priori_IEG_{subclass}.csv"))
  
  # continue if nothing significant; nothing to plot
  if (nrow(res |> filter(padj<0.05)) == 0) {
    next
  }
  
  # gene types ----------------------------------------
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
      linewidth = 1 
    ) +
    facet_wrap(vars(gene), ncol=3, scales = 'free_y') +
    labs(title = glue('ZT:EE genes in {subclass}'),
         x = 'Pseudobulk samples',
         y = 'Expression') +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          strip.background = element_blank(),
          legend.position = 'none')
  print(p1)
  
  
  # line plot ----------------------------------------
  p2 <- df_means |> 
    ggplot() +
    aes(x = ZT, y = mean_count, color = activity_condition, group = activity_condition) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = mean_count - sem, ymax = mean_count + sem), width = 0.1) +
    facet_wrap(vars(gene), ncol=3, scales = "free_y") +
    scale_color_manual(values = activity_colors) +
    labs(title = glue('ZT:EE genes in {subclass}'),
         y = 'Expression') +
    theme_classic() +
    theme(axis.text.x = element_text(size = 8, angle = 30),
          axis.title.x = element_blank(),
          axis.line.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title = element_text(size = 15),
          plot.title = element_text(size = 20, hjust = 0.5),
          strip.background = element_blank(),
          legend.position = 'none')
  print(p2)
  
  
  # save plots ----------------------------------------
  if (SAVE_PLOTS) {
    subclass_fname <- ShrinkSubclassName(subclass)
    save_dir = '05-results/ORCA/raw_R_plots'
    
    save_dims <- tibble(
      celltype = c('CA1', 'DG'),
      width = c(5, 3),
      height = c(4, 1.5)
    ) |> 
      filter(celltype == subclass) |> 
      select(-celltype) |> 
      as.list() |> 
      unlist()
    
    # png
    ggsave(plot = p1,
           path = save_dir,
           filename = glue('genomewide__{subclass_fname}_barplot.png'),
           width =  5, height = 4, dpi = 900, bg = 'white')
    ggsave(plot = p2,
           path = save_dir,
           filename = glue('genomewide__{subclass_fname}_lineplot.png'),
           width =  3.4, height = 1.5, dpi = 900, bg = 'white')
    # svg
    ggsave(plot = p1 + LoadBarebonesTheme(ticks = 'y'),
           path = save_dir,
           filename = glue('genomewide_{subclass_fname}_barplot.svg'),
           width =  6, height = 4)
    ggsave(plot = p2 + LoadBarebonesTheme(ticks = 'both'),
           path = save_dir,
           filename = glue('genomewide_{subclass_fname}_lineplot.svg'),
           width =  6, height = 4)
  }
}
