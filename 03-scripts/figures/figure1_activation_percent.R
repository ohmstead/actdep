# ---- load libs & data ----
print("Loading libraries and data...")

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(patchwork)
library(glue)

source("03-scripts/R/seq_functions.R")
nuclei <- LoadDataset("May2024", "combined")
activity_colors <- LoadActivityColors("Dec2024")


# ---- get active cells IEG ----
print("Getting percent of cells active using IEGs...")

IEG_symbols <- LoadGeneList("IEG")

outputs <- FindActiveCells(nuclei, gene_list = IEG_symbols, gene_threshold = 3)
df_active_cells_IEGs <- outputs$df_active_cells
activation_thresholds_IEGs <- outputs$activation_thresholds


# ---- plot % active IEGs ----
print("Plotting percent of cells active using IEGs...")

p_activation_IEG <- ggplot(df_active_cells_IEGs, aes(x = activity_condition, y = num_upregd_genes, fill = activity_condition)) +
  geom_jitter(width = 0.3, height = 0.25, shape = 21, alpha = 0.3) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 2.5, linewidth = 0.5, linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = activity_colors) +
  theme_minimal() +
  # make y-axis labels at 0, 3, 5, and 10
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

lncRNA_degs <- read_csv("04-analysis/DEGs/May2024_activity_condition/016_CA1-ProS_Glut_EE30m_vs_SE.csv") |> 
  left_join(gene_biotypes, by = c('gene' = 'name')) |> 
  filter(biotype == 'lncRNA') |> 
  filter(chrom != 'mm39_X' & chrom != 'mm39_Y') |> 
  filter(abs(avg_log2FC) > 0.5)

outputs <- FindActiveCells(nuclei, gene_list = lncRNA_degs$gene, gene_threshold = 2)
df_active_cells_lncRNA <- outputs$df_active_cells
activation_thresholds_lncRNA <- outputs$activation_thresholds


# ---- plot % active lncRNA ----
print("Plotting percent of cells active using lncRNAs...")

p_activation_lncRNA <- ggplot(df_active_cells_lncRNA, aes(x = activity_condition, y = num_upregd_genes, fill = activity_condition)) +
  geom_jitter(width = 0.3, height = 0.25, shape = 21, size = 1, alpha = 0.5) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 1.5, linewidth = 0.5, linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = activity_colors) +
  theme_minimal() +
  scale_y_continuous(breaks = c(0, 2, 5, 10)) +
  coord_cartesian(ylim = c(0, 6)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        legend.position = 'none') +
  labs(x = '', y = '', fill = 'Condition')
print(p_activation_lncRNA)
print(lncRNA_degs$gene)


# ---- plot KDE inset ----
print("Plotting KDE example insets...")

p_SE <- df_active_cells_IEGs |> 
  filter(activity_condition == 'SE') |> 
ggplot() +
  geom_density(aes(x = Npas4, fill = activity_condition)) +
  geom_vline(xintercept = 0.5, color = 'black', linetype = 'dotted', lwd=1) +
  geom_label(x = 0.55, y = 2.5, label = '90th expression percentile in HC cells', size = 6, color = 'black', hjust = 0) +
  geom_label(x = 0.52, y = 1, label = "'active'", size = 4, color = "black", hjust = 0) +
  geom_label(x = 0.48, y = 1, label = "'inactive'", size = 4, color = "black", hjust = 1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  labs(title = "Npas4 in HC cells", x = "Npas4 expression") +
  theme_bw() +
  theme(
    legend.position = 'none',
  )

p_EE_30m <- df_active_cells_IEGs |> 
  filter(activity_condition == 'EE30m') |> 
  ggplot() +
  geom_density(aes(x = Npas4, fill = activity_condition)) +
  geom_vline(xintercept = 0.5, color = 'black', linetype = 'dotted', lwd=1) +
  geom_label(x = 0.55, y = 1.5, label = '90th expression percentile in HC cells', size = 6, color = 'black', hjust = 0) +
  geom_label(x = 0.52, y = 0.75, label = "'active'", size = 4, color = "black", hjust = 0) +
  geom_label(x = 0.48, y = 0.75, label = "'inactive'", size = 4, color = "black", hjust = 1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  labs(title = "Npas4 in KA cells", x = "Npas4 expression") +
  theme_bw() +
  theme(
    legend.position = 'none',
  )

p_KA <- df_active_cells_IEGs |> 
  filter(activity_condition == 'KA') |> 
  ggplot() +
  geom_density(aes(x = Npas4, fill = activity_condition)) +
  geom_vline(xintercept = 0.5, color = 'black', linetype = 'dotted', lwd=1) +
  geom_label(x = 0.55, y = 0.65, label = '90th expression percentile in HC cells', size = 6, color = 'black', hjust = 0) +
  geom_label(x = 0.52, y = 0.25, label = "'active'", size = 4, color = "black", hjust = 0) +
  geom_label(x = 0.48, y = 0.25, label = "'inactive'", size = 4, color = "black", hjust = 1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  labs(title = "Npas4 in 30m EE cells", x = "Npas4 expression") +
  theme_bw() +
  theme(
    legend.position = 'none',
  )

print(p_SE / p_EE_30m / p_KA)


# ---- save plots ----
if (SAVE_PLOTS) {
  print("Saving plots...")
  ggsave(plot = p_activation_IEG,
         path = '05-results/figure1/raw_R_plots/', 
         filename = 'activation_IEGs.png', 
         width = 7, height = 3, dpi = 900)
  ggsave(plot = p_activation_lncRNA,
         path = '05-results/figure1/raw_R_plots/', 
         filename = 'activation_lncRNAs.png', 
         width = 7, height = 3, dpi = 900)
  # inset plots
  ggsave(plot = p_SE, 
         path = '05-results/figure1/raw_R_plots/', 
         filename = 'HC_inset.png', 
         width = 8, height = 2.5, dpi = 900)
  ggsave(plot = p_EE_30m,
         path = '05-results/figure1/raw_R_plots/', 
         filename = 'EE_30m_inset.png', 
         width = 8, height = 2, dpi = 900)
  ggsave(plot = p_KA,
         path = '05-results/figure1/raw_R_plots/', 
         filename = 'KA_inset.png', 
         width = 8, height = 2.5, dpi = 900)
} else {
  print("Plotted without saving...")
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))