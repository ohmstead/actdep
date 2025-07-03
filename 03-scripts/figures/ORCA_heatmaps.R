# Makes hm plots for Clock & IEG expression in SE subclasses
# Makes hm plots for significant interaction genes in CA1 and DG

library(DESeq2)
library(ComplexHeatmap)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset("Dec2024")
zt_colors <- LoadZTColors()
activity_colors <- LoadActivityColors()

# define genes to test ----------------------------------------
a_priori_genes <- c(LoadGeneList('circadian'), LoadGeneList('IEG'))

subclass_list <- LoadSubclassesToUse(nuclei)

seurat_subsets <- list(
  CA1   = nuclei |> subset(subclass_name == subclass_list[1]  & activity_condition %in% c('SE', 'EE30m')),
  DG    = nuclei |> subset(subclass_name == subclass_list[7]  & activity_condition %in% c('SE', 'EE30m')),
  Astro = nuclei |> subset(subclass_name == subclass_list[17] & activity_condition %in% c('SE', 'EE30m')),
  Oligo = nuclei |> subset(subclass_name == subclass_list[18] & activity_condition %in% c('SE', 'EE30m'))
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
for (subclass in names(dds_subsets)) {
  # get hm count data
  dds <- dds_subsets[[subclass]]
  
  bulk_counts <- counts(dds, normalized = T) |> 
    t() |> 
    as.data.frame() |>
    rownames_to_column('sample') |> 
    select(sample, all_of(LoadGeneList('IEG'))) |> 
    tibble() |> 
    print()
  
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
    mutate(sample = factor(sample, levels = unique(sample))) |> 
    print()
  
  # 2) find average FC of EE30m from SE
  df_means <- df_bulk |>
    group_by(gene, activity_condition, ZT) |>
    summarize(mean_count = mean(count, na.rm = T)+1, .groups = "drop") |> 
    group_by(gene, ZT) |> 
    summarize(log2FC = log2(
      mean_count[activity_condition == 'EE30m'] / 
      mean_count[activity_condition == 'SE']
      )) |>
    print()

  # 3) plot
  p <- df_means |>
    mutate(gene = factor(gene, levels = rev(LoadGeneList('IEG')))) |> 
  ggplot() +
    aes(x = ZT, y = gene, fill = log2FC) +
    geom_tile(color = 'black', linewidth = 0.25) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                         midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
    labs(
      title = paste0(subclass),
      x = 'ZT',
      y = 'Gene',
      fill = 'log2FC'
    ) +
    theme(axis.title = element_blank(),
          axis.text.x  = element_text(size = 20, angle = 30, vjust = 0.8),
          axis.text.y  = element_text(size = 20, face = 'italic'),
          plot.title = element_text(size = 25, hjust = 0.5),
          panel.grid.major = element_blank(),
          # legend.position = 'none'
          )
  print(p)
}


# Clock hm ----------------------------------------
genes_to_plot <- rev(c('Per1', 'Per2', 'Cry2', 'Clock', 'Bmal1'))
for (subclass in names(dds_subsets)) {
  # get hm count data
  dds <- dds_subsets[[subclass]]
  
  bulk_counts <- counts(dds, normalized = T) |> 
    t() |> 
    as.data.frame() |>
    rownames_to_column('sample') |> 
    select(sample, all_of(genes_to_plot)) |> 
    tibble() |> 
    print()
  
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
    mutate(sample = factor(sample, levels = unique(sample))) |> 
    print()
  
  # 2) find average FC of EE30m from SE
  df_means <- df_bulk |>
    group_by(gene, activity_condition, ZT) |>
    summarize(mean_count = mean(count, na.rm = T)+1, .groups = "drop") |> 
    group_by(gene, ZT) |> 
    summarize(log2FC = log2(
      mean_count[activity_condition == 'EE30m'] / 
        mean_count[activity_condition == 'SE']
    )) |>
    mutate(gene = factor(gene, levels = genes_to_plot)) |> 
    print()
  
  # 3) plot
  p <- ggplot(df_means) +
    aes(x = ZT, y = gene, fill = log2FC) +
    geom_tile(color = 'black', linewidth = 0.25) +
    scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                         midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
    labs(
      title = paste0(subclass),
      x = 'ZT',
      y = 'Gene',
      fill = 'log2FC'
    ) +
    theme(axis.title = element_blank(),
          axis.text.x  = element_text(size = 20, angle = 30, vjust = 0.8),
          axis.text.y  = element_text(size = 20, face = 'italic'),
          plot.title = element_text(size = 25, hjust = 0.5),
          panel.grid.major = element_blank(),
          legend.position = 'none'
    )
  print(p)
}


# tmp work with Brenda ----------------------------------------
subclass <- 'DG'
subclass <- 'CA1'
dds <- dds_subsets[[subclass]]

bulk_counts <- counts(dds, normalized = T) |> 
  t() |> 
  as.data.frame() |>
  rownames_to_column('sample') |> 
  # select(sample, all_of(c('Cecr2', 'Daglb', 'Enoph1', 'Fcho2', 'Gm5820', 'Itpkb', 'Poc1b', 'Trmt61a', 'Zfp420', 'Zfp523'))) |> 
  select(sample, all_of(c('Kcnn2', 'Napepld'))) |>
  tibble() |> 
  print()

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
  mutate(sample = factor(sample, levels = unique(sample))) |> 
  print()

# 2) find average FC of EE30m from SE
df_means <- df_bulk |>
  group_by(gene, activity_condition, ZT) |>
  summarize(mean_count = mean(count, na.rm = T)+1, .groups = "drop") |> 
  group_by(gene, ZT) |> 
  summarize(log2FC = log2(
    mean_count[activity_condition == 'EE30m'] / 
      mean_count[activity_condition == 'SE']
  )) |>
  print()

p <- ggplot(df_means) +
  aes(x = ZT, y = gene, fill = log2FC) +
  geom_tile(color = 'black', linewidth = 0.25) +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', 
                       midpoint = 0, limits = c(-3, 3), oob = scales::squish) +
  labs(
    title = paste0(subclass),
    x = 'ZT',
    y = 'Gene',
    fill = 'log2FC'
  ) +
  theme(axis.title = element_blank(),
        axis.text.x  = element_text(size = 20, angle = 30, vjust = 0.8),
        axis.text.y  = element_text(size = 20, face = 'italic'),
        plot.title = element_text(size = 25, hjust = 0.5),
        panel.grid.major = element_blank(),
        legend.position = 'none'
  )
print(p)
