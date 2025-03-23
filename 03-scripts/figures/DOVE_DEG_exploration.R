# This script finds DEGs at 30m and 6h for astrocytes. 
# It characterizes upregulated biological processes using GO analysis.
library(ggplot2)
library(patchwork)
library(glue)
library(readr)
library(dplyr)

library(Seurat)

source('03-scripts/R/seq_functions.R')

activity_colors <- LoadActivityColors()
supertype_colors <- LoadAllenColors('supertype')


# read in DEGs
read_DEG_file <- function(fname) {
  read_csv(fname, show_col_types = F) |> 
    # filter(padj < 0.05, abs(log2FoldChange) > 0.585) |> 
    arrange(desc(log2FoldChange)) |> 
    mutate(contrast = str_sub(fname,-15,-5),
           de = ifelse((padj<0.05&abs(log2FoldChange)>0.585), T, F)) |> 
    print()
}
deg_EE30m <- read_DEG_file("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/319_Astro-TE_NN__EE30m_vs_SE.csv")
deg_EE6h  <- read_DEG_file("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/319_Astro-TE_NN__EE6h_vs_SE.csv")
deg_KA30m <- read_DEG_file("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/319_Astro-TE_NN__KA30m_vs_SE.csv")
deg_KA6h  <- read_DEG_file("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/319_Astro-TE_NN__KA6h_vs_SE.csv")

deg <- tibble(rbind(deg_EE30m, deg_EE6h, deg_KA30m, deg_KA6h)) |> print()

# get bg genes expressed in at least 5% of astrocytes
astro <- nuclei |> subset(subclass_name == '319 Astro-TE NN')
mat <- astro[["SCT"]]@data
mat <- mat[rowMeans(mat > 0) > 0.01, ]
bg_genes <- rownames(mat)

# get GO terms
ego_EE30m <- RunGOEnrichment(deg_EE30m |> filter(de==T) |> pull(gene), bg_genes)
ego_EE6h  <- RunGOEnrichment(deg_EE_6h |> filter(de==T) |> pull(gene), bg_genes)
ego_KA30m <- RunGOEnrichment(deg_KA_30m |> filter(de==T) |> pull(gene), bg_genes)
ego_KA6h  <- RunGOEnrichment(deg_KA_6h |> filter(de==T) |> pull(gene), bg_genes)

p <- ggplot(deg_KA6h) +
  aes(log2FoldChange, -log10(padj), color = de) +
  geom_point(aes(color = de)) +
  geom_hline(yintercept = 0.05, linetype = 'dashed')
HoverLocator(p, information = deg_KA6h)


# specific genes ---------------------------------------------
go_1904893 <- c(
  'Cav1',
  'Ggnbp2',
  'Sac3d1',
  'Hnf4a',
  'Hmga2',
  'Sh2b3',
  'Cish',
  'Neurod1',
  'Tbx1',
  'Bcl3',
  'Nf2',
  'Inpp5f',
  'Ptprd',
  'Ptprc',
  'Leprot',
  'Suz12',
  'Parp14',
  'Socs1',
  'Dab1',
  'Adipor1',
  'Hgs',
  'Pibf1',
  'Socs3',
  'Socs2',
  'Traf3ip1',
  'Ptpn2',
  'Gbp7',
  'Irf1'
)
VlnPlot(astro, features = go_1904893, split.by = 'activity_condition', ncol = 10)
# GO:1904893
# GO:0046426

genes_to_plot <- deg_EE30m |> 
  filter(de == T) |>
  pull(gene)

for (gene in genes_to_plot) {
  p1 <- FeaturePlot(astro, gene, pt.size = 1, alpha=1) + scale_color_viridis_c(option = 'inferno')
  p2 <- VlnPlot(astro, gene, split.by = 'activity_condition')
  p <- p1 + p2
  ggsave(glue('04-analysis/astrocytes/{gene}.png'), plot=p, width = 7, height = 10)
}


# look at IMN distribution over activity conditions
imn <- nuclei |> subset(subclass_name == '038 DG-PIR Ex IMN')
imn@meta.data |> 
ggplot() +
  aes(x = activity_condition, fill = supertype_name) +
  geom_bar(position = 'stack') +
  scale_fill_manual(values = supertype_colors)


# just make a damn dotplot ---------------------------------------------
interesting_feats <- c(
  '9330159F19Rik', 'Ddah1', 'Dio2', 'Mxi1', 'Pde10a', 'Epas1', 'Plekho2',
  'Fam20a', 'Slco1c1', 'Foxo1', 'Usp2', 'Hdac4', 'Zswim6', 'Kcnn2', 'Map3k19'
)
astro |> 
  RidgePlot(features = interesting_feats, group.by = 'activity_condition', cols = activity_colors)
  scale_fill_manual(values = activity_colors)
astro |> 
  DotPlot(features = interesting_feats, group.by = 'activity_condition', cols = 'viridis')

FeaturePlot(astro, interesting_feats) +
  scale_color_viridis_c()
