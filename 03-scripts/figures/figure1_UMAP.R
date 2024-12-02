# ---- load libs & data ----
print("Loading libraries and data...")
library(Seurat)
library(ggplot2)

source('03-scripts/R/seq_functions.R')
nuclei <- LoadDataset("May2024")

condition_colors <- LoadConditionColors("May2024")
subclass_colors <- LoadAllenColors(clade = "subclass")

# ---- plot ----
print("Plotting UMAPs...")

p1 <- DimPlot(nuclei, reduction = "umap", group.by = "subclass_name", pt.size = 2, shuffle = TRUE, seed = 17) +
  scale_color_manual(values = subclass_colors) +
  theme_minimal() +
  labs(title = '', x = '', y = '') +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = 'none'
  )

p2 <- DimPlot(nuclei, reduction = "umap", group.by = "condition", pt.size = 2, shuffle = TRUE, seed = 17) +
  scale_color_manual(values = condition_colors) +
  theme_minimal() +
  labs(title = '', x = '', y = '') +
  theme(
    panel.grid = element_blank(),
    axis.text = element_blank(),
    legend.text = element_text(size = 12),
    legend.position = 'none'
  )

p1
p2

# ---- save plot ----
if (SAVE_PLOTS) {
  print("Saving plot...")
  
  ggsave(path = '05-results/figure1/raw_R_plots/', filename = 'UMAP_subclass.png', plot = p1, width = 3, height = 3, dpi = 900)
  ggsave(path = '05-results/figure1/raw_R_plots/', filename = 'UMAP_condition.png', plot = p2, width = 3, height = 3, dpi = 900)
} else {
  print("Plotting without saving...")
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))