## ----Fig4A-B
# This script plots the percent of active cells in each condition for the CA1 subclass.
# It does no analysis at the supertype level. It uses the canonical IEG list.

library(ggridges)

source("03-scripts/R/seq_functions.R")

if (!exists('nuclei')) {nuclei <- LoadDataset("Dec2024")}

activity_colors <- LoadActivityColors("Dec2024", palette = 'grayscale')


# get active cells  ----------------------------------------
print("Getting percent of cells active using IEGs...")

IEG_symbols <- LoadGeneList("IEG")
if (!exists('FindActiveCells_CA1_outputs')) {
  FindActiveCells_CA1_outputs <- FindActiveCells(
    nuclei, 
    subclass = '016 CA1-ProS Glut', 
    gene_list = IEG_symbols, 
    gene_threshold = 3
  )
}
df_active_cells_IEGs <- FindActiveCells_CA1_outputs$df_active_cells


# plot activation  ----------------------------------------
p_activation_IEG <- ggplot(df_active_cells_IEGs) +
  aes(x = activity_condition, y = num_upregd_genes, fill = activity_condition) +
  geom_jitter(alpha = 0.3, size = 0.25, width = 0.3, height = 0.25, shape = 21, set.seed(17)) +
  geom_boxplot(width = 0.3, outlier.shape = NA) +
  geom_hline(yintercept = 2.5, linewidth = 0.5, linetype = 'dashed', color = 'black') +
  scale_fill_manual(values = activity_colors) +
  scale_y_continuous(breaks = c(0, 3, 5, 10)) +
  labs(x = '', y = 'Num. IEGs upregulated', title = 'CA1-ProS activation') +
  theme(
    panel.grid = element_blank(),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 14)
  )

# print(p_activation_IEG)


# example insets  ----------------------------------------
p_inset_Fos <- df_active_cells_IEGs |> 
  filter(activity_condition %in% c('SE', 'EE30m', 'KA30m')) |>
  mutate(activity_condition = factor(activity_condition, levels = rev(c('SE', 'EE30m', 'KA30m')))) |> 
ggplot() +
  aes(x = Fos, y = activity_condition, fill = activity_condition) +
  geom_density_ridges(scale = 3, alpha = 0.8) +
  geom_vline(xintercept = FindActiveCells_CA1_outputs$activation_thresholds['Fos'], color = 'black', linetype = 'dotted', lwd=1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 3), expand = FALSE) +
  labs(title = 'Fos') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 12, face='italic'),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )

p_inset_Arc <- df_active_cells_IEGs |> 
  filter(activity_condition %in% c('SE', 'EE30m', 'KA30m')) |>
  mutate(activity_condition = factor(activity_condition, levels = rev(c('SE', 'EE30m', 'KA30m')))) |> 
ggplot() +
  aes(x = Arc, y = activity_condition, fill = activity_condition) +
  geom_density_ridges(scale = 3, alpha = 0.8) +
  geom_vline(xintercept = FindActiveCells_CA1_outputs$activation_thresholds['Arc'], color = 'black', linetype = 'dotted', lwd=1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 3), expand = FALSE) +
  labs(title = 'Arc') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 12, face='italic'),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )

p_inset_Nr4a1 <- df_active_cells_IEGs |> 
  filter(activity_condition %in% c('SE', 'EE30m', 'KA30m')) |>
  mutate(activity_condition = factor(activity_condition, levels = rev(c('SE', 'EE30m', 'KA30m')))) |> 
ggplot() +
  aes(x = Nr4a1, y = activity_condition, fill = activity_condition) +
  geom_density_ridges(scale = 3, alpha = 0.8) +
  geom_vline(xintercept = FindActiveCells_CA1_outputs$activation_thresholds['Nr4a1'], color = 'black', linetype = 'dotted', lwd=1) +
  scale_fill_manual(values = activity_colors) +
  coord_cartesian(xlim = c(0, 3.5), expand = FALSE) +
  labs(title = 'Nr4a1') +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 12, face='italic'),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank()
  )

layout <- c("
AD
BD
CD
")
p_inset_Fos + p_inset_Arc + p_inset_Nr4a1 + p_activation_IEG + 
  plot_layout(design=layout, widths = c(1,3))


# save plots ----------------------------------------
if (SAVE_PLOTS) {
  print("Saving plots...")
  save_path <- "05-results/Figure4/raw_R_plots"
  # PNG
  ggsave(plot = p_activation_IEG,
         path = save_path, 
         filename = 'activation_IEGs.png', 
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
## ----