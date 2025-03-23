library(ggplot2)
library(dplyr)
library(glue)
library(tibble)
library(rstatix)

library(Seurat)

source('03-scripts/R/seq_functions.R')
save_path <- "05-results/ORCA/raw_R_plots"

activity_colors <- LoadActivityColors('Dec2024')
# nuclei <- LoadDataset('Dec2024')
# nuclei_subclass <- nuclei |> subset(subclass_name == '016 CA1-ProS Glut' & activity_condition == 'SE')


calc_sem <- function(data, varname, groupnames) {
  #+++++++++++++++++++++++++
  # Function to calculate the mean and the standard deviation
  # for each group
  #+++++++++++++++++++++++++
  # data : a data frame
  # varname : the name of a column containing the variable
  #to be summariezed
  # groupnames : vector of column names to be used as
  # grouping variables
  require(plyr)
  
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sem = sd(x[[col]], na.rm=TRUE) / sqrt(length(x[[col]])) )
  }
  
  data_sum <- ddply(data, groupnames, .fun=summary_func, varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# SE ZT insets ----------------------------------------
known_circadian_genes <- LoadGeneList('circadian')
mat_circadian <- GetAssayData(nuclei_subclass)
mat_circadian <- mat_circadian[known_circadian_genes,] |> t()
nuclei_subclass@meta.data <- nuclei_subclass@meta.data |> 
  mutate(ZT.collection = case_when(
    (ZT == 'ZT0' & activity_condition == 'EE6h') ~ 'ZT6',
    (ZT == 'ZT4' & activity_condition == 'EE6h') ~ 'ZT10',
    (ZT == 'ZT12' & activity_condition == 'EE6h') ~ 'ZT18',
    (ZT == 'ZT16' & activity_condition == 'EE6h') ~ 'ZT22',
    TRUE ~ ZT
  )) |> 
  mutate(ZT.collection = factor(ZT.collection, 
                                levels = c('ZT0', 'ZT4', 'ZT6', 'ZT10', 'ZT12', 'ZT16', 'ZT18', 'ZT22')))

# merge circadian gene expression with cell metadata
df_circadian <- as.matrix(mat_circadian) |> 
  as.data.frame() |> 
  rownames_to_column(var = 'barcode') |> 
  as_tibble() |> 
  right_join(nuclei_subclass@meta.data) |> 
  tidyr::pivot_longer(cols = all_of(known_circadian_genes), names_to = 'gene', values_to = 'expression')


for (gene_to_plot in known_circadian_genes) {
  p <- df_circadian |> 
    calc_sem('expression', c('ZT', 'gene')) |>
    filter(gene == gene_to_plot) |> 
  ggplot() +
    aes(x = ZT, y = expression, group = gene) +
    geom_line(linewidth = 2) +
    geom_errorbar(aes(ymin = expression-sem, ymax = expression+sem), linewidth = 2, width = 0.1) +
    geom_point(size = 5) +
    labs(title = glue('{gene_to_plot} in CA1 SE'),
         y = 'Normalized expression') +
    theme(
      axis.text.x = element_text(size = 15, angle = 30),
      axis.title.y = element_text(size = 15),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 15, hjust = 0.5)
    )
  print(p)
  
  if (SAVE_PLOTS) {
    save_path <- "05-results/MOTH/raw_R_plots"
    ggsave(plot = p,  # png
           filename = glue('CA1_{gene_to_plot}_inset.png'),
           path = save_path,
           width = 4, height = 5, units = 'in', dpi = 900, bg = 'white')
    ggsave(plot = p + LoadBarebonesTheme(ticks = 'y'),  # svg
           filename = glue('CA1_{gene_to_plot}_inset.svg'),
           path = save_path,
           width = 5, height = 5, units = 'in')
  }
}


# VlnPlots ----------------------------------------
known_circadian_genes <- LoadGeneList('circadian')

ca1 <- subset(nuclei, subclass_name == '016 CA1-ProS Glut')
mat_circadian <- GetAssayData(ca1)
mat_circadian <- mat_circadian[known_circadian_genes,] |> t()
df_circadian <- as.matrix(mat_circadian) |> 
  as.data.frame() |> 
  rownames_to_column(var = 'barcode') |> 
  as_tibble() |> 
  right_join(ca1@meta.data) |> 
  tidyr::pivot_longer(cols = all_of(known_circadian_genes), names_to = 'gene', values_to = 'expression')

for (gene_to_plot in known_circadian_genes) {
  # Compute maximum expression per ZT for the gene and conditions of interest.
  max_expr <- df_circadian %>%
    filter(gene == gene_to_plot, activity_condition %in% c("SE", "EE30m")) %>%
    group_by(ZT) %>%
    summarise(y_max = max(expression, na.rm = TRUE))
  
  # Compute the t-test results per ZT and join with the maximum expression.
  stat_table <- df_circadian %>% 
    filter(gene == gene_to_plot, activity_condition %in% c("SE", "EE30m")) %>%
    group_by(ZT) %>%
    wilcox_test(expression ~ activity_condition) %>% 
    adjust_pvalue(method = "bonferroni") %>% 
    add_significance("p.adj") %>% 
    mutate(x = as.numeric(factor(ZT)),
           xmin = x - 0.15,
           xmax = x + 0.15) %>%
    left_join(max_expr, by = "ZT") %>%
    # Set y.position to be 5% above the maximum expression in that ZT group.
    mutate(y.position = y_max * 1.1)
  
  p <- df_circadian |> 
    filter(gene == gene_to_plot, activity_condition %in% c('SE', 'EE30m')) |> 
    ggplot() +
    aes(x = ZT, y = expression, fill = activity_condition) +
    geom_point(shape = 21, size = 5, alpha = 0.4, 
               position = position_jitterdodge(dodge.width = 0.75, 
                                               jitter.width = 1, 
                                               jitter.height = 0.05,
                                               seed = 17)
    ) +
    geom_violin(aes(group = interaction(ZT, activity_condition)),
                fill = 'gray75', width = 0.5,
                position = position_dodge(width = 0.75)) +
    stat_pvalue_manual(stat_table, 
                       label = "p.adj.signif", 
                       inherit.aes = FALSE) +
    scale_fill_manual(values = activity_colors) +
    theme_minimal() +
    theme(
      legend.position = 'none',
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 20),
      panel.border = element_rect(linewidth = 1, fill = NA)
    )

  print(p)
  
  if (SAVE_PLOTS) {
    ggsave(plot = p,  # png
           filename = glue('CA1_{gene_to_plot}_VlnPlot.png'),
           path = save_path,
           width = 4, height = 5, units = 'in', dpi = 900, bg = 'white')
    ggsave(plot = p + LoadBarebonesTheme(ticks = 'y'),  # svg
           filename = glue('CA1_{gene_to_plot}_VlnPlot.svg'),
           path = save_path,
           width = 4, height = 5, units = 'in')
  }
}


# DotPlots ----------------------------------------
for (gene_to_plot in known_circadian_genes) {
  # DotPlot
  p <- df_circadian |> 
    filter(gene == gene_to_plot, activity_condition %in% c('SE', 'EE30m')) |> 
    group_by(activity_condition, ZT) |> 
    dplyr::summarize(
      pct.expressed = mean(expression > 0),
      expression = mean(expression)
    ) |>
    ggplot() +
    aes(y = activity_condition, x = ZT, color = expression, size = pct.expressed) +
    geom_point(shape = 20) +
    scale_size_continuous(
      name = 'Percent cells expressing',
      range = c(1, 10),
      breaks = c(0.1, 0.25, 0.5, 0.75, 1),
      guide = 'legend'
    ) +
    scale_color_viridis_c(option = 'D', limits = c(0, 2)) +
    theme_minimal() +
    labs(title = glue('{gene_to_plot}')) +
    theme(
      legend.position = 'bottom',
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      axis.text = element_text(size = 20),
      panel.border = element_rect(linewidth = 1, fill = NA)
    )
  print(p)
  
  if (SAVE_PLOTS) {
    ggsave(plot = p,  # png
           filename = glue('CA1_{gene_to_plot}_DotPlot.png'),
           path = save_path,
           width = 4, height = 5, units = 'in', dpi = 900, bg = 'white')
    ggsave(plot = p + LoadBarebonesTheme(ticks = 'y'),  # svg
           filename = glue('CA1_{gene_to_plot}_DotPlot.svg'),
           path = save_path,
           width = 4, height = 5, units = 'in')
  }
}


# facet plots (subclass) ----------------------------------------
known_circadian_genes <- c(
  'Per1', 'Per2', 'Per3', 'Clock', 'Bmal1', 'Cry1', 'Cry2',
  'Nr1d1', 'Nr1d2', 'Hif3a', 'Thbs3', 'Hspa5', 'Fkbp5'
)

subclasses <- c(
  '016 CA1-ProS Glut',
  # '017 CA3 Glut',
  '037 DG Glut',
  '319 Astro-TE NN',
  # '326 OPC NN',
  '327 Oligo NN'
  # '334 Microglia NN'
)

for (subclass in subclasses) {
  nuclei_tmp <- nuclei |>
    subset(subclass_name == subclass & activity_condition %in% c('SE'))

  mat_circadian <- GetAssayData(nuclei_tmp)
  mat_circadian <- mat_circadian[LoadGeneList("IEG"),] |> t()
  nuclei_tmp@meta.data <- nuclei_tmp@meta.data |>
    mutate(ZT.collection = case_when(
      (ZT == 'ZT0' & activity_condition == 'EE6h') ~ 'ZT6',
      (ZT == 'ZT4' & activity_condition == 'EE6h') ~ 'ZT10',
      (ZT == 'ZT12' & activity_condition == 'EE6h') ~ 'ZT18',
      (ZT == 'ZT16' & activity_condition == 'EE6h') ~ 'ZT22',
      TRUE ~ ZT
    )) |>
    mutate(ZT.collection = factor(ZT.collection,
                                  levels = c('ZT0', 'ZT4', 'ZT6', 'ZT10', 'ZT12', 'ZT16', 'ZT18', 'ZT22')))

  # merge circadian gene expression with cell metadata
  df_circadian <- as.matrix(mat_circadian) |>
    as.data.frame() |>
    rownames_to_column(var = 'barcode') |>
    as_tibble() |>
    right_join(nuclei_tmp@meta.data) |>
    tidyr::pivot_longer(cols = all_of(LoadGeneList("IEG")), names_to = 'gene', values_to = 'expression')

  p <- df_circadian |>
    calc_sem('expression', c('ZT', 'gene')) |>
  ggplot() +
    aes(x = ZT, y = expression, group = gene) +
    geom_line(linewidth = 2) +
    geom_errorbar(aes(ymin = expression-sem, ymax = expression+sem), linewidth = 2, width = 0.1) +
    geom_point() +
    facet_wrap(~gene, scales = 'free_y') +
    labs(title = subclass) +
    theme(
      strip.text = element_text(size = 30),
      axis.title = element_blank(),
      axis.text = element_text(size = 20),
      plot.title = element_text(size = 30)
    )
  print(p)
}
