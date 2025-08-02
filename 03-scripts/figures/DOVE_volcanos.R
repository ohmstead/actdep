source('03-scripts/R/seq_functions.R')

library(ggrepel)

activity_colors  <- LoadActivityColors()
volcano_colors <- list(
  upregulated = '#E85D75',
  downregulated = '#5AA9E6',
  no_change   = 'gray80'
)

# read in DEGs ----------------------------------------
degs_30m <- list(
  Astro = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/319_Astro-TE_NN__EE30m_vs_SE.csv"),
  Oligo = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/327_Oligo_NN__EE30m_vs_SE.csv"),
  OPC   = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/326_OPC_NN__EE30m_vs_SE.csv"),
  Micro = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/334_Microglia_NN__EE30m_vs_SE.csv")
)

degs_6h <- list(
  Astro = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/319_Astro-TE_NN__EE6h_vs_SE.csv"),
  Oligo = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/327_Oligo_NN__EE6h_vs_SE.csv"),
  OPC   = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/326_OPC_NN__EE6h_vs_SE.csv"),
  Micro = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/334_Microglia_NN__EE6h_vs_SE.csv")
)


# plot volcanos ----------------------------------------
plots_30m <- list()
plots_6h  <- list()

for (celltype in names(degs_30m)) {
  deg_30m <- degs_30m[[celltype]]
  deg_6h  <- degs_6h[[celltype]]
  
  p1 <- ggplot(deg_30m) +
    aes(x = log2FoldChange.shrink, y = -log10(padj), color = classification) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = 'dashed', color = 'black') +
    geom_text_repel(aes(label = ifelse(classification != 'no_change', gene, '')), size = 3) +
    scale_color_manual(values = volcano_colors) +
    labs(title = glue("{celltype} EE30m vs SE")) +
    theme_minimal() +
    theme(legend.position = "none")
  p1
  plots_30m[[celltype]] <- p1
  
  p2 <- ggplot(deg_6h) +
    aes(x = log2FoldChange.shrink, y = -log10(padj), color = classification) +
    geom_point() +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed') +
    geom_vline(xintercept = c(-0.585, 0.585), linetype = 'dashed', color = 'black') +
    geom_text_repel(aes(label = ifelse(classification != 'no_change', gene, '')), size = 3) +
    scale_color_manual(values = volcano_colors) +
    labs(title = glue("{celltype} EE6h vs SE")) +
    theme_minimal() +
    theme(legend.position = "none")
  p2
  plots_6h[[celltype]] <- p2
  
  print(p1 + p2)
}


# save plots ----------------------------------------
save_path <- "05-results/DOVE/raw_R_plots/"

if (SAVE_PLOTS) {
  # 30m volcanos
  for (p in plots_30m) {
    # png
    ggsave(plot = p,
           filename = glue("volcano__30m_{p}.png"),
           path = save_path,
           width = 4, height = 3)
    # svg
    ggsave(plot = p + LoadBarebonesTheme(ticks = 'both'),
           filename = glue("volcano_30m_{p}.svg"),
           path = save_path,
           width = 4, height = 3)
  }
  
  # 6h volcanos
  for (p in plots_6h) {
    # png
    ggsave(plot = p,
           filename = glue("volcano__6h_{p}.png"),
           path = save_path,
           width = 4, height = 3)
    # svg
    ggsave(plot = p + LoadBarebonesTheme(ticks = 'both'),
           filename = glue("volcano_6h_{p}.svg"),
           path = save_path,
           width = 4, height = 3)
  }
}