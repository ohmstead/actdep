# ---- load libs & data ----
print("Loading libraries and data...")

library(Seurat)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(patchwork)
library(glue)

source("03-scripts/R/seq_functions.R")
nuclei <- LoadDataset("Dec2024")
activity_colors <- LoadActivityColors("Dec2024", palette = 'grayscale')


# ---- get active cells ----
print("Getting percent of cells active using IEGs...")

IEG_symbols <- LoadGeneList("IEG")
lncRNA_symbols <- LoadGeneList("lncRNA")

outputs_IEG    <- FindActiveCells(nuclei, subclass = '016 CA1-ProS Glut', gene_list = IEG_symbols, gene_threshold = 3)
outputs_lncRNA <- FindActiveCells(nuclei, subclass = '016 CA1-ProS Glut', gene_list = lncRNA_symbols, gene_threshold = 3)

df_active_cells_IEGs <- outputs_IEG$df_active_cells
df_active_cells_lncRNA <- outputs_lncRNA$df_active_cells


# ---- plot % active IEGs ----
print("Plotting percent of cells active using IEGs...")

p_activation_IEG <- ggplot(df_active_cells_IEGs) +
  aes(x = activity_condition, y = num_upregd_genes, fill = activity_condition) +
  geom_jitter(alpha = 0.3, size = 0.25, width = 0.3, height = 0.25, shape = 21, set.seed(17)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_hline(yintercept = 2.5, linewidth = 0.5, linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = activity_colors) +
  theme_minimal() +
  scale_y_continuous(breaks = c(0, 3, 5, 10)) +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none'
  ) +
  labs(x = '', y = '', fill = 'Condition')

print(p_activation_IEG)
print(IEG_symbols)


# ---- plot % active lncRNA ----
print("Plotting percent of cells active using lncRNA...")

p_activation_lncRNA <- ggplot(df_active_cells_lncRNA) +
  aes(x = activity_condition, y = num_upregd_genes, fill = activity_condition) +
  geom_jitter(alpha = 0.3, size = 0.25, width = 0.3, height = 0.25, shape = 21, set.seed(17)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
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

print(p_activation_lncRNA)
print(lncRNA_symbols)


# ---- KDE insets ----
print("Plotting KDE example insets...")

p_inset_Fos <- df_active_cells_IEGs |> 
  filter(activity_condition %in% c('SE', 'EE30m', 'KA30m')) |>
  mutate(activity_condition = factor(activity_condition, levels = rev(c('SE', 'EE30m', 'KA30m')))) |> 
ggplot() +
  aes(x = Fos, y = activity_condition, fill = activity_condition) +
  geom_density_ridges(scale = 3, alpha = 0.8) +
  geom_vline(xintercept = outputs_IEG$activation_thresholds['Fos'], color = 'black', linetype = 'dotted', lwd=1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  theme_void() +
  theme(legend.position = 'none')
print(p_inset_Fos)

p_inset_Arc <- df_active_cells_IEGs |> 
  filter(activity_condition %in% c('SE', 'EE30m', 'KA30m')) |>
  mutate(activity_condition = factor(activity_condition, levels = rev(c('SE', 'EE30m', 'KA30m')))) |> 
ggplot() +
  aes(x = Arc, y = activity_condition, fill = activity_condition) +
  geom_density_ridges(scale = 3, alpha = 0.8) +
  geom_vline(xintercept = outputs_IEG$activation_thresholds['Arc'], color = 'black', linetype = 'dotted', lwd=1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  theme_void() +
  theme(legend.position = 'none')
print(p_inset_Arc)

p_inset_Nr4a1 <- df_active_cells_IEGs |> 
  filter(activity_condition %in% c('SE', 'EE30m', 'KA30m')) |>
  mutate(activity_condition = factor(activity_condition, levels = rev(c('SE', 'EE30m', 'KA30m')))) |> 
ggplot() +
  aes(x = Nr4a1, y = activity_condition, fill = activity_condition) +
  geom_density_ridges(scale = 3, alpha = 0.8) +
  geom_vline(xintercept = outputs_IEG$activation_thresholds['Nr4a1'], color = 'black', linetype = 'dotted', lwd=1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 5), expand = FALSE) +
  theme_void() +
  theme(legend.position = 'none')
print(p_inset_Nr4a1)


# ---- save plots ----
if (SAVE_PLOTS) {
  print("Saving plots...")
  save_path <- "05-results/LION/raw_R_plots"
  # PNG
  ggsave(plot = p_activation_IEG,
         path = save_path, 
         filename = 'activation_IEGs.png', 
         width = 7, height = 3, dpi = 300)
  ggsave(plot = p_activation_lncRNA,
         path = save_path, 
         filename = 'activation_lncRNAs.png', 
         width = 7, height = 3, dpi = 300)
  # inset plots
  ggsave(plot = p_inset_Fos, 
         path = save_path, 
         filename = 'inset_Fos.png', 
         width = 8, height = 2.5, dpi = 900)
  ggsave(plot = p_inset_Arc, 
         path = save_path, 
         filename = 'inset_Arc.png', 
         width = 8, height = 2.5, dpi = 900)
  ggsave(plot = p_inset_Nr4a1, 
         path = save_path, 
         filename = 'inset_Nr4a1.png', 
         width = 8, height = 2.5, dpi = 900)
  # SVG
  ggsave(plot = p_activation_IEG,
         path = save_path, 
         filename = 'activation_IEGs.svg', 
         width = 7, height = 3)
  ggsave(plot = p_activation_lncRNA,
         path = save_path, 
         filename = 'activation_lncRNAs.svg', 
         width = 7, height = 3)
  ggsave(plot = p_inset_Fos, 
         path = save_path, 
         filename = 'inset_Fos.svg', 
         width = 8, height = 2.5)
  ggsave(plot = p_inset_Arc, 
         path = save_path, 
         filename = 'inset_Arc.svg', 
         width = 8, height = 2.5)
  ggsave(plot = p_inset_Nr4a1, 
         path = save_path, 
         filename = 'inset_Nr4a1.svg', 
         width = 8, height = 2.5)
} else {
  print("Plotted without saving...")
}
