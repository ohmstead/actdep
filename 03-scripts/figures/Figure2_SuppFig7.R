## ----IMN volcano: SE vs EE30m and SE vs EE6h
source('03-scripts/R/seq_functions.R')

library(ggrepel)

# load DEG csvs ----------------------------------------
csv_paths <- list(
  EE30m = "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/038_DG-PIR_Ex_IMN__EE30m_vs_SE.csv",
  EE6h = "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/038_DG-PIR_Ex_IMN__EE6h_vs_SE.csv"
)
csv_results <- lapply(csv_paths, read_csv)

# define colormap for points
cmap <- list('downregulated' = hcl.colors(n=5, palette = 'RdBu')[4],
             'no_change' = '#C2B8B2',
             'upregulated' = hcl.colors(n=5, palette = 'RdBu')[2])

# EE30m vs SE ----
p_EE30m <- csv_results[['EE30m']] |> ggplot() +
  aes(x = log2FoldChange.shrink,
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
  geom_text_repel(data = csv_results[['EE30m']] |>
                    filter(classification != 'no_change'),
                  inherit.aes = T,
                  max.overlaps = 10) +
  scale_fill_manual(values = cmap) +
  scale_color_manual(values = cmap) +
  labs(
    title = paste0('IMN: EE30m vs SE'),
    x = 'log2(FC from SE)',
    y = '-log10(FDR)') +
  scale_y_continuous(limits = c(0, -log10(0.001))) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15),
    # remove minor gridlines
    panel.grid.minor = element_blank(),
  )


# EE6h vs SE ----
p_EE6h <- csv_results[['EE6h']] |> ggplot() +
  aes(x = log2FoldChange.shrink,
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
  geom_text_repel(data = csv_results[['EE6h']] |>
                    filter(classification != 'no_change'),
                  inherit.aes = T,
                  max.overlaps = 10) +
  scale_fill_manual(values = cmap) +
  scale_color_manual(values = cmap) +
  labs(
    title = paste0('IMN: EE6h vs SE'),
    x = 'log2(FC from SE)',
    y = '') +
  scale_y_continuous(limits = c(0, -log10(0.001))) +
  scale_x_continuous(limits = c(-1.5,1.5)) +
  theme(
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.title = element_text(size = 15),
    # remove minor gridlines
    panel.grid.minor = element_blank(),
  )


p <- Reduce(`+`, list(p_EE30m, p_EE6h)) + plot_layout(nrow = 1)
print(p)


# save volcanos ----
if (exists('SAVE_PLOTS') & SAVE_PLOTS) {
  save_path = "05-results/Figure2_SuppFig3/"
  ggsave(
    filename = "Volcano_IMN_vs_SE.png",
    plot = p,
    path = save_path,
    width = 12, height = 4, dpi = 300
  )
  # svgsave(
  #   filename = "Volcano_IMN_vs_SE.svg",
  #   plot = p + LoadBarebonesTheme(ticks = 'both'),
  #   path = save_path,
  #   width = 12, height = 4
  # )
}


# GO analysis SE vs EE6h ----
# target genes = significant DEGs; background = all genes tested (non-NA padj)
df_EE6h <- csv_results$EE6h

target_genes     <- df_EE6h |> filter(classification != 'no_change') |> pull(gene)
background_genes <- df_EE6h |> filter(!is.na(padj)) |> pull(gene)

GO_results <- RunGOEnrichment(
  target_gene_symbols = target_genes,
  background_gene_list = background_genes
) |> arrange(pvalue)

# NOTE: RunGOEnrichment() returns clusterProfiler's full tested-term table
# (its `@result` slot is unfiltered), not just terms passing the pvalueCutoff/
# qvalueCutoff used internally (p.adjust < 0.05, qvalue < 0.2) — so check
# significance explicitly rather than relying on nrow(). With only 12 target
# genes (8 up, 4 down), no GO BP terms actually survive that FDR cutoff here;
# fall back to showing the top terms by nominal p-value, flagged as such below.
is_FDR_significant <- any(GO_results$p.adjust < 0.05 & GO_results$qvalue < 0.2)