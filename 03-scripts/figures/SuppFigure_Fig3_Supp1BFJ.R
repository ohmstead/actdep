## ----SuppFig3_Supp1BFJ
source('03-scripts/R/seq_functions.R')

library(ggrepel)

# load DEG csvs ----------------------------------------
csv_paths <- list(
  CA1 = "04-analysis/reviewer_response/voom/016_CA1-ProS_Glut__EE30m_vs_SE.csv",
  SST = "04-analysis/reviewer_response/voom/053_Sst_Gaba__EE30m_vs_SE.csv",
  Astrocyte = "04-analysis/reviewer_response/voom/319_Astro-TE_NN__EE30m_vs_SE.csv"
)
csv_results <- lapply(csv_paths, read_csv)

# make plots ----------------------------------------
# define colormap for points
cmap <- list('downregulated' = hcl.colors(n=5, palette = 'RdBu')[4],
             'no_change' = '#C2B8B2',
             'upregulated' = hcl.colors(n=5, palette = 'RdBu')[2])

plots <- list()

for (celltype in names(csv_results)) {
  df <- csv_results[[celltype]]

  # voom's unshrunk logFC estimates are not bounded the way apeglm shrinkage
  # bounds them, so compute data-driven symmetric x-limits that never clip
  # a data point off the panel (while still spanning at least [-1.5, 4],
  # matching the original DESeq2-based plot's range)
  lim_val <- max(abs(range(c(-1.5, 4, df$log2FoldChange), na.rm = TRUE)))

  p <- df |> ggplot() +
    aes(x = log2FoldChange,
        y = -log10(padj),
        fill = classification,
        color = classification,
        label = gene) +
    geom_point(
      shape  = 21,
      size   = 3,
      color = 'black',
      stroke = 0.1
    ) +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = 'dashed') +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_text_repel(data = df |>
                      filter(classification != 'no_change'),
                    inherit.aes = T,
                    max.overlaps = 10) +
    scale_fill_manual(values = cmap) +
    scale_color_manual(values = cmap) +
    labs(
      title = paste0(celltype, ': EE30m vs SE'),
      x = 'log2(FC from SE)',
      y = '-log10(FDR)') +
    scale_x_continuous(limits = c(-lim_val, lim_val), breaks = pretty(c(-lim_val, lim_val))) +
    scale_y_continuous(breaks = seq(0, ceiling(max(-log10(df$padj), na.rm = TRUE)), by = 5)) +
    theme(
      legend.position = 'none',
      plot.title = element_text(hjust = 0.5, size = 15),
      axis.title = element_text(size = 15),
      # remove minor gridlines
      panel.grid.minor = element_blank(),
    )
  # print(p)

  plots[[celltype]] <- p
}

p <- Reduce(`+`, plots) + plot_layout(nrow = 1)
print(p)


# save ----------------------------------------
if (exists('SAVE_PLOTS') & SAVE_PLOTS) {
  save_path = "05-results/SuppFigure_Fig3_Supp1/raw_R_plots"
  for (celltype in names(plots)) {
    # png
    ggsave(
      filename = glue("Volcano__{celltype}_EE30m_vs_SE.png"),
      plot = plots[[celltype]],
      path = save_path,
      width = 6, height = 4, dpi = 300
    )

    # svg
    svgsave(
      filename = glue("Volcano_{celltype}_EE30m_vs_SE.svg"),
      plot = plots[[celltype]] + LoadBarebonesTheme(ticks = 'both'),
      path = save_path,
      width = 9, height = 3
    )
  }
}
## ----
