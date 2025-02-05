# ---- load libs & data ----
print("Loading libraries and data...")

library(Seurat)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(glue)

source("03-scripts/R/seq_functions.R")
nuclei <- LoadDataset("Dec2024")
activity_colors <- LoadActivityColors("Dec2024")


# ---- get active cells IEG ----
print("Getting percent of cells active using IEGs...")

IEG_symbols <- LoadGeneList("IEG")

outputs <- FindActiveCells(nuclei, gene_list = IEG_symbols, gene_threshold = 3)
df_active_cells_IEGs <- outputs$df_active_cells
activation_thresholds_IEGs <- outputs$activation_thresholds


# ---- plot % active IEGs ----
print("Plotting percent of cells active using IEGs...")

p_activation_IEG <- ggplot(df_active_cells_IEGs) +
  aes(x = activity_condition, y = num_upregd_genes, fill = activity_condition) +
  geom_jitter(width = 0.3, height = 0.25, shape = 21, alpha = 0.3, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 2.5, linewidth = 0.5, linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = activity_colors) +
  theme_minimal() +
  scale_y_continuous(breaks = c(0, 3, 5, 10)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'none'
  ) +
  labs(x = '', y = '', fill = 'Condition')

print(p_activation_IEG)
print(IEG_symbols)


# ---- get active cells lncRNA ----
print("Getting percent of cells active using lncRNAs...")

gene_biotypes <- read_csv("04-analysis/all_gene_stats.csv")

lncRNA_list <- LoadGeneList('lncRNA')

lnc_gene_thresh <- round(0.2 * length(lncRNA_list$gene))
# lnc_gene_thresh <- 4
outputs <- FindActiveCells(nuclei, gene_list = lncRNA_list$gene, gene_threshold = lnc_gene_thresh)
df_active_cells_lncRNA <- outputs$df_active_cells
activation_thresholds_lncRNA <- outputs$activation_thresholds


# ---- plot % active lncRNA ----
print("Plotting percent of cells active using lncRNAs...")

p_activation_lncRNA <- ggplot(df_active_cells_lncRNA) +
  aes(x = activity_condition, y = num_upregd_genes, fill = activity_condition) +
  geom_jitter(width = 0.3, height = 0.25, shape = 21, alpha = 0.3, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 3-0.5, linewidth = 0.5, linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = activity_colors) +
  theme_minimal() +
  scale_y_continuous(breaks = c(0, 3, 5, 10)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'none'
  ) +
  labs(x = '', y = '', fill = 'Condition')

print(p_activation_lncRNA)
print(lncRNA_degs$gene)


# ---- plot KDE inset ----
print("Plotting KDE example insets...")

p_inset <- df_active_cells_IEGs |> 
  filter(activity_condition %in% c('SE', 'EE30m', 'KA30m')) |>
  mutate(activity_condition = factor(activity_condition, levels = rev(c('SE', 'EE30m', 'KA30m')))) |> 
ggplot() +
  aes(x = Npas4, y = activity_condition, fill = activity_condition) +
  geom_density_ridges(scale = 3, alpha = 0.8) +
  geom_vline(xintercept = 0.5, color = 'black', linetype = 'dotted', lwd=1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  theme_void() +
  theme(legend.position = 'none')

print(p_inset)


# ---- save plots ----
if (SAVE_PLOTS) {
  print("Saving plots...")
  ggsave(plot = p_activation_IEG,
         path = '05-results/figure2/raw_R_plots/', 
         filename = 'activation_IEGs.png', 
         width = 7, height = 3, dpi = 900)
  ggsave(plot = p_activation_lncRNA,
         path = '05-results/figure2/raw_R_plots/', 
         filename = 'activation_lncRNAs.png', 
         width = 7, height = 3, dpi = 900)
  # inset plot
  ggsave(plot = p_inset, 
         path = '05-results/figure2/raw_R_plots/', 
         filename = 'inset.png', 
         width = 8, height = 2.5, dpi = 900)
} else {
  print("Plotted without saving...")
}
