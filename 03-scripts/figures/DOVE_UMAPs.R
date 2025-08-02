source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset('Dec2024')

activity_colors <- LoadActivityColors()
subclass_colors <- LoadAllenColors()


# make UMAPs ----------------------------------------
DimPlot(nuclei, 
        cells.highlight = WhichCells(nuclei, idents = c('319 Astro-TE NN')),
        cols.highlight  = subclass_colors['319 Astro-TE NN'],
        pt.size = 0.3,
        sizes.highlight = 1) +
  theme_void() +
  theme(legend.position = 'none')

# save
ggsave(plot = p_UMAPs,
       filename = glue("UMAPs_319_Astro-TE_NN.svg"),
       path = save_path,
       width = 4, height = 4)