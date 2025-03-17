# ---- load libs & data ----
print("Loading libraries and data...")

library(Seurat)
library(tidyverse)
library(patchwork)
library(glue)
library(ggh4x)

source("03-scripts/R/seq_functions.R")

nuclei <- LoadDataset("Dec2024")

activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors('subclass')
subclass_list <- LoadSubclassesToUse(nuclei, ascertainment = 'custom')
df_degs <- read_csv('04-analysis/DEGs/Dec2024_activity_condition_minPct5/0_sorted_DEG_list.csv')

# get active cells in each subclass ----
print("Getting percent of cells active using IEGs...")

df_percent_active_all = tibble()
df_threshold_finding_all = tibble()
counter = 1

for(current_subclass in subclass_list) {
  print(glue("{current_subclass} ({counter}/{length(subclass_list)})"))

  # get genes specific to this subclass
  subclass_genes <- df_degs |> 
    filter(n_subclasses == 1) |> 
    filter(subclass == current_subclass) |> 
    filter(classification == 'ERG') |> 
    pull(gene)

  # find the threshold based on the number of genes in the subclass-specific list
  outputs <- FindActiveCells(nuclei, subclass = current_subclass, gene_list = subclass_genes)
  
  # extract and aggregate outputs
  df_active_cells   <- outputs$df_active_cells

  # define range of thresholds to test
  threshold_values <- seq(1, max(df_active_cells$num_upregd_genes, na.rm = TRUE), by = 1)

  # calculate percent activated at each threshold
  df_threshold_finding <- map_df(threshold_values, function(thresh) {
    percent_active <- df_active_cells |> 
      group_by(activity_condition) |> 
      summarize(percent_active = mean(num_upregd_genes >= thresh)) |> 
      mutate(threshold = thresh) |> 
      mutate(subclass = current_subclass)
  })

  # find the threshold where the percent active is closest to target_pct
  target_pct <- 0.15
  target_threshold <- df_threshold_finding |> 
    filter(activity_condition == 'SE') |> 
    filter(abs(percent_active - target_pct) == min(abs(percent_active - target_pct))) |> 
    pull(threshold)

  df_percent_active_final <- df_threshold_finding |> 
    filter(threshold == target_threshold) |> 
    print()
  
  df_threshold_finding_all <- bind_rows(df_threshold_finding_all, df_threshold_finding)
  df_percent_active_all    <- bind_rows(df_percent_active_all, df_percent_active_final)

  counter <- counter+1
}

df_percent_active_all$subclass <- factor(df_percent_active_all$subclass, levels = subclass_list)









# tests with Brenda ----
pv <- subset(nuclei, subclass_name == '052 Pvalb Gaba')
sst <- subset(nuclei, subclass_name == '053 Sst Gaba')
cck <- subset(nuclei, subclass_name == '047 Sncg Gaba')
VlnPlot(pv, 'Pvalb', group.by = 'activity_condition')
VlnPlot(sst, 'Sst', group.by = 'activity_condition')
VlnPlot(cck, 'Cck', group.by = 'activity_condition')

# tests with Danny ----
# re-run DEG for microglia
microglia <- subset(nuclei, subclass_name == '334 Microglia NN')
Idents(microglia) <- microglia$activity_condition

microglia_degs.05 <- FindMarkers(
  microglia, 
  logfc.threshold = 0.585,
  test.use = 'MAST',
  min.pct = 0.05,
  ident.1 = 'EE30m', ident.2 = 'SE',
  latent.vars = c('sublibrary', 'percent.mt'),
) |> 
  rownames_to_column(var = 'gene') |> 
  filter(p_val_adj < 0.05) |> 
  arrange(desc(avg_log2FC)) |> 
  print()

microglia_thresholds <- tibble('gene' = names(outputs$activation_thresholds), 
                               'expression_threshold' = outputs$activation_thresholds)

tmp <- df_active_cells |>
  filter(activity_condition == 'EE30m') |>
  pivot_longer(cols = -c(cell, activity_condition, active_binary, num_upregd_genes), names_to = 'gene', values_to = 'expression') |> 
  left_join(microglia_thresholds, by = 'gene') |> 
  mutate(binary_expression = expression > 0) |> 
  group_by(gene) |> 
  mutate(pct_expressing = mean(binary_expression)) |> 
  summarize(pct_expressing = first(pct_expressing)) |> 
  arrange(desc(pct_expressing)) |> 
  mutate(gene = factor(gene, levels = gene))


ggplot(tmp) +
  aes(y = gene, x = pct_expressing) +
  geom_col() +
  geom_vline(xintercept = 0.005, color = 'red') +
  geom_hline(yintercept = 45)

# explore the expression distribution for a highly expressed gene
df_active_cells |> 
  filter(!str_detect(activity_condition, 'KA')) |>
ggplot() +
  aes(x = Agpat4, fill = activity_condition) +
  geom_density(alpha = 0.3) +
  geom_histogram(aes(y = ..count../sum(..count..)), position = 'dodge') +
  scale_fill_manual(values = activity_colors) +
  scale_x_continuous(limits = c(0, 3))


# se_df <- df_percent_active_all |>
  # filter(activity_condition == "SE") |>
  # select(threshold, subclass, se_percent = percent_active)

# Join back to the original data frame and compute the fold-change.
df_threshold_finding_all |>
  group_by(threshold, subclass) |> 
  mutate(fold_change = percent_active / first(percent_active[activity_condition == "SE"])) |> 
  mutate(threshold = factor(threshold)) |> 
filter(activity_condition == 'EE30m') |> 
ggplot() +
  aes(x = threshold, y = fold_change, fill = threshold) +
  geom_col(position = 'dodge') +
  facet_wrap(~subclass)
  scale_fill_brewer()





# plot activation within each subclass ----
subclass_colors_strip <- subclass_colors[names(subclass_colors) %in% subclass_list]
subclass_colors_strip <- subclass_colors_strip[match(subclass_list, names(subclass_colors_strip))]
# strip = strip_themed(background_x = elem_list_rect(fill = subclass_colors_strip))

p_subclass <- ggplot(df_percent_active_all) +
  aes(x = activity_condition, y = percent_active, fill = activity_condition) +
  geom_col(position = 'dodge') +
  geom_text(aes(label = scales::percent(percent_active, accuracy = 1),
                color = activity_condition,
                size = 8), 
            vjust = -0.2,  # Positions text slightly above bars
            size = 4) +
  # facet_wrap2(~subclass, strip = strip, nrow = 4, ncol = 5) +
  facet_wrap2(~subclass, nrow = 4, ncol = 5) +
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = activity_colors) +
  scale_color_manual(values = activity_colors) +
  theme(legend.position = 'none',
        strip.text = element_text(size = 12),
        axis.title.x = element_blank(),
        axis.title.y = element_blank())


# save plots ----
if (SAVE_PLOTS) {
  print("Saving plots...")
  ggsave(plot = p_subclass,
         path = '05-results/figure2/raw_R_plots/', 
         filename = 'activation_by_subclass.png', 
         width = 18, height = 12, dpi = 900)
} else {
  print("Plotted without saving...")
}


# examine dentate cells bc they're fishy
outputs <- FindActiveCells(nuclei, subclass = '037 DG Glut', gene_list = IEG_symbols, gene_threshold = 3)
dentate <- outputs$df_active_cells

ggplot(dentate) +
  aes(y = num_upregd_genes, fill = activity_condition) +
  geom_histogram(position = 'dodge')

mat <- GetAssayData(nuclei, 'SCT', layer = 'data')

corr_mat <- GetCorrData(nuclei, 'SE', IEG_symbols)
hm <- PlotComplexHeatmap(corr_mat, plot_title = 'Dentate EE 30m')
draw(hm)

mat <- mat[
  rownames(mat) %in% IEG_symbols,  # get genes
  colnames(mat) %in% dentate$cell   # get dentate cells
] |> 
as.matrix()

hm <- ComplexHeatmap::Heatmap(t(mat))

# get SE and EE30m cell names
cells.30m <- dentate |> 
  filter(activity_condition == 'EE30m') |> 
  pull(cell)
cells.se <- dentate |> 
  filter(activity_condition == 'SE') |> 
  pull(cell)

# plot for SE and 30m independently
mat.30m <- mat[,colnames(mat) %in% cells.30m]
mat.se <- mat[,colnames(mat) %in% cells.se]

hm.30m <- ComplexHeatmap::Heatmap(t(mat.30m))
hm.se  <- ComplexHeatmap::Heatmap(t(mat.se))
draw(hm.30m)
draw(hm.se)