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
IEG_symbols <- LoadGeneList("IEG")

# get active cells in each subclass  ----------------------------------------
print("Getting percent of cells active using IEGs...")

df_percent_active_all = tibble()
counter = 1

for(subclass in subclass_list) {
  print(glue("{subclass} ({counter}/{length(subclass_list)})"))

  outputs <- FindActiveCells(nuclei, subclass = subclass, gene_list = IEG_symbols, gene_threshold = 3)
  
  # extract and aggregate outputs
  df_active_cells   <- outputs$df_active_cells
  df_percent_active <- outputs$df_percent_active
  df_percent_active_all <- bind_rows(df_percent_active_all, df_percent_active)

  counter <- counter+1
}

df_percent_active_all$subclass <- factor(df_percent_active_all$subclass, levels = subclass_list)


# plot activation within each subclass  ----------------------------------------
subclass_colors_strip <- subclass_colors[names(subclass_colors) %in% subclass_list]
subclass_colors_strip <- subclass_colors_strip[match(subclass_list, names(subclass_colors_strip))]
strip = strip_themed(background_x = elem_list_rect(fill = subclass_colors_strip))

p <- ggplot(df_percent_active_all) +
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


# save plots  ----------------------------------------
if (SAVE_PLOTS) {
  print("Saving plots...")
  save_path <- "05-results/LION/raw_R_plots"
  ggsave(plot = p,
         path = save_path, 
         filename = 'activation_subclass_IEG.png', 
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