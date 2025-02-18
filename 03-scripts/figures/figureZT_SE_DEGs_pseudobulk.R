# Load required libraries
library(Seurat)
library(DESeq2)
library(ggrepel)
library(tidyverse)
library(SingleCellExperiment)
library(ComplexHeatmap)
source("03-scripts/R/seq_functions.R")

# define fxns ----
getConstrastResults <- function(dds, contrast, threshold) {
  res <- results(dds, contrast = contrast) |> 
    as.data.frame() |> 
    rownames_to_column("gene") |> 
    relocate(gene) |> 
    mutate(
      classification = 
        ifelse(log2FoldChange > thresh & padj < 0.05, "upregulated", 
               ifelse(log2FoldChange < -thresh & padj < 0.05, "downregulated", 
                      "no_change")
        )
    )
    
    print(tibble(res))

  rownames(res) <- res$gene
  
  return(tibble(res))
}

volcanoPlot <- function(res, title, threshold) {
  point_color <- setNames(c("red", "blue", "black"), 
                          c("upregulated", "downregulated", "no_change"))
  
  p <- ggplot(res) +
    aes(x = log2FoldChange, y = -log10(padj), color = classification) +
    geom_point() +
    geom_text_repel(
      data = subset(res, classification != "no_change"),  # Only label DEGs
      aes(label = gene),
      size = 3,
      max.overlaps = 20,
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-threshold, threshold), linetype = "dashed") +
    scale_color_manual(values = point_color) +
    labs(x = "log2FoldChange", y = "-log10(padj)", title = title) +
    theme(legend.position = "none")
  
  print(p)
  return(p)
}


# load data ----
nuclei <- LoadDataset('Dec2024')
celltype <- subset(nuclei, 
                   (subclass_name == '016 CA1-ProS Glut' &
                    activity_condition == 'SE'))

# prep data ----
# scrambled_indices <- sample(seq_len(nrow(celltype@meta.data)))   # scramble cell identities!!
# celltype@meta.data$sample <- celltype@meta.data$sample[scrambled_indices]
# celltype@meta.data$ZT <- celltype@meta.data$ZT[scrambled_indices]

# only test genes expressed in 5% of cells
mat <- celltype[["SCT"]]@data
mat <- mat[rowMeans(mat > 0) > 0.05, ]
cell_metadata <- celltype@meta.data
gene_list <- rownames(mat)

celltype$Condition <- factor(celltype$activity_condition, levels = c("SE", "EE30m", "EE6h"))
celltype$ZT <- factor(celltype$ZT, levels = c("ZT0", "ZT4", "ZT12", "ZT16"))
celltype$sample <- factor(celltype$sample)



pseudobulk_counts <- AggregateExpression(celltype, group.by = "sample", features = gene_list, return.seurat = FALSE)$RNA

col_data <- celltype@meta.data |> 
  as_tibble() |> 
  distinct(sample, Condition, ZT) |> 
  mutate(sample = str_replace(sample, '_', '-')) |> 
  column_to_rownames("sample") |> 
  filter(!is.na(Condition)) |> 
  arrange(ZT, Condition)

col_data <- col_data[colnames(pseudobulk_counts),]
col_data$sample <- rownames(col_data)


# run DE ----
dds <- DESeqDataSetFromMatrix(countData = as.matrix(pseudobulk_counts), 
                              colData = col_data, 
                              design = ~ ZT)
dds <- DESeq(dds)

thresh = 0.585

# ZT4 vs ZT0 ----
contrast = c("ZT", "ZT4", "ZT0")
res_ZT4_vs_ZT0 <- getConstrastResults(dds, contrast, thresh)
p <- volcanoPlot(res_ZT4_vs_ZT0, "ZT4 vs ZT0", thresh)


# ZT12 vs ZT0 ----
contrast = c("ZT", "ZT12", "ZT0")
res_ZT12_vs_ZT0 <- getConstrastResults(dds, contrast, thresh)
p <- volcanoPlot(res_ZT12_vs_ZT0, "ZT12 vs ZT0", thresh)


# ZT16_vs_ZT0 ----
contrast = c("ZT", "ZT16", "ZT0")
res_ZT16_vs_ZT0 <- getConstrastResults(dds, contrast, thresh)
p <- volcanoPlot(res_ZT16_vs_ZT0, "ZT16 vs ZT0", thresh)


# ZT12 vs ZT4 ----
contrast = c("ZT", "ZT12", "ZT4")
res_ZT12_vs_ZT4 <- getConstrastResults(dds, contrast, thresh)
p <- volcanoPlot(res_ZT12_vs_ZT4, "ZT12 vs ZT4", thresh)


# ZT16 vs ZT4 ----
contrast = c("ZT", "ZT16", "ZT4")
res_ZT16_vs_ZT4 <- getConstrastResults(dds, contrast, thresh)
volcanoPlot(res_ZT16_vs_ZT4, "ZT16 vs ZT4", thresh)


# ZT16 vs ZT12 ----
contrast = c("ZT", "ZT16", "ZT12")
res_ZT16_vs_ZT12 <- getConstrastResults(dds, contrast, thresh)
p <- volcanoPlot(res_ZT16_vs_ZT12, "ZT16 vs ZT12", thresh)


# heatmap ----
# assemble all genes and their contrast of origin
all_genes <- bind_rows(
  res_ZT4_vs_ZT0 |> mutate(contrast = "ZT4 vs ZT0") |> filter(classification != 'no_change'),
  res_ZT12_vs_ZT0 |> mutate(contrast = "ZT12 vs ZT0") |> filter(classification != 'no_change'),
  res_ZT16_vs_ZT0 |> mutate(contrast = "ZT16 vs ZT0") |> filter(classification != 'no_change'),
  res_ZT12_vs_ZT4 |> mutate(contrast = "ZT12 vs ZT4") |> filter(classification != 'no_change'),
  res_ZT16_vs_ZT4 |> mutate(contrast = "ZT16 vs ZT4") |> filter(classification != 'no_change'),
  res_ZT16_vs_ZT12 |> mutate(contrast = "ZT16 vs ZT12") |> filter(classification != 'no_change')
)

# remove duplicates
gene_list <- all_genes |> 
  distinct(gene) |> 
  pull()

# get DESeq-normalized expression
mat <- counts(dds, normalized = T)
gene_mat <- mat[gene_list,] |> 
  t() |> 
  as.data.frame() |> 
  rownames_to_column('sample') |> 
  pivot_longer(cols = -sample, names_to = 'gene', values_to = 'expression') |> 
  mutate(ZT = str_extract(sample, "ZT\\d+")) |> 
  mutate(ZT = factor(ZT, levels = c("ZT0", "ZT4", "ZT12", "ZT16"))) |> 
  group_by(gene, ZT) |> 
  mutate(mean_expression = mean(expression)) |> 
  ungroup() |> 
  group_by(gene) |> 
  mutate(z_expression = scale(mean_expression)) |> 
  select(gene, ZT, z_expression) |> 
  distinct() |> 
  pivot_wider(names_from = ZT, values_from = z_expression) |> 
  column_to_rownames('gene') |> 
  as.matrix()

k = 4
col_fun <- circlize::colorRamp2(c(-2,0,2), hcl_palette = "blue-red2")
hm <- Heatmap(gene_mat,
              col = col_fun,
              cluster_columns = FALSE,
              row_split = k,
              cluster_row_slices = F
              # row_names_gp = gpar(fontsize = 5)
              )
hm = draw(hm)

InteractiveComplexHeatmap::htShiny(hm)


# find p-values for all gene lists  ----


if (SAVE_PLOTS) {
  png("05-results/figureZT/raw_R_plots/SE_ZT_DEGs.png", width = 5, height = 10, units = "in", res = 900)
  print(hm)
  dev.off()
}

# GO analysis ----
cluster_genes <- lapply(row_order(hm), function(idx) rownames(mat)[idx])

go_terms_all <- tibble()

terms.1 <- RunGOEnrichment(cluster_genes[[1]]) |> mutate(cluster = 1)
go_terms_all <- rbind(go_terms_all, terms.1)

terms.2 <- RunGOEnrichment(cluster_genes[[2]]) |> mutate(cluster = 2)
go_terms_all <- rbind(go_terms_all, terms.2)

terms.3 <- RunGOEnrichment(cluster_genes[[3]]) |> mutate(cluster = 3)
go_terms_all <- rbind(go_terms_all, terms.3)

terms.4 <- RunGOEnrichment(cluster_genes[[4]]) |> mutate(cluster = 4)
go_terms_all <- rbind(go_terms_all, terms.4)

terms.5 <- RunGOEnrichment(cluster_genes[[5]]) |> mutate(cluster = 5)
go_terms_all <- rbind(go_terms_all, terms.5)