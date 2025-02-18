# Load required libraries
library(Seurat)
library(DESeq2)
library(readr)
library(ggrepel)
library(tidyverse)
library(MAST)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(rpca)
source("03-scripts/R/seq_functions.R")

ZT_colors <- LoadZTColors()


# define fxns ----
volcanoPlot <- function(res, title) {
  point_color <- setNames(c("red", "blue", "black"), c("upregulated", "downregulated", "no_change"))
  
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
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    scale_color_manual(values = point_color) +
    labs(x = "log2FoldChange", y = "-log10(padj)", title = title) +
    theme(legend.position = "none")
  
  print(p)
  return(p)
}

plotPCA <- function(sca_obj){
  set.seed(17)
  projection <- rpca(t(assay(sca_obj)))$x
  colnames(projection)=c("PC1","PC2","PC3","PC4")
  pca <- data.table(projection,  as.data.frame(colData(sca_obj)))
  print(ggpairs(pca, columns=c('PC1', 'PC2', 'PC3', 'libSize', 'PercentToHuman', 'nGeneOn', 'exonRate'),
          mapping=aes(color=condition), upper=list(continuous='blank')))
  invisible(pca)
}

quickViolin <- function(genes) {
  p <- VlnPlot(celltype, group.by = 'ZT', features = genes) +
    scale_fill_manual(values = ZT_colors)
  print(p)

  return(p)
}


# load data ----
valid_conditions <- c('SE', 'EE30m', 'EE6h')
nuclei <- LoadDataset('Dec2024')
celltype <- subset(
  nuclei, 
  (
    subclass_name == '016 CA1-ProS Glut' & 
    # activity_condition %in% valid_conditions
    activity_condition == 'SE'
  )
)

# convert to SingleCellAssay
mat <- celltype[["SCT"]]@data

# remove rows from mat that have nonzeros in fewer than 5% of cells
mat <- mat[rowMeans(mat > 0) > 0.05, ]
cell_metadata <- celltype@meta.data
gene_metadata <- data.frame(gene = rownames(mat))
sca <- FromMatrix(exprsArray = as.matrix(mat),
                  cData = cell_metadata,
                  fData = gene_metadata)

# calcualte the cellular detection rate after gene filt
cdr <-colSums(assay(sca)>0)
colData(sca)$cdr <- scale(cdr)

# save genes tested
genes_tested <- rownames(mat)
csv_write(genes_tested, "04-analysis/DEGs/Dec2024_ZT/genes_tested.csv")


# viz data ----
# plotPCA(sca)


# run DE with ZT0 as the base level ----
options(mc.cores = 16)
zlm_model.ZT0 <- zlm(~ZT + cdr + (1|sample), 
                 data = sca,
                 method = 'glmer',
                 ebayes = FALSE,
                 parallel = TRUE)


# run DE with ZT4 as the base level ----
sca$ZT <- factor(sca$ZT, levels = c("ZT4", "ZT0", "ZT12", "ZT16"))
zlm_model.ZT4 <- zlm(~ZT + cdr + (1|sample), 
                 data = sca,
                 method = 'glmer',
                 ebayes = FALSE,
                 parallel = TRUE)


# run DE with ZT12 as the base level ----
sca$ZT <- factor(sca$ZT, levels = c("ZT12", "ZT0", "ZT4", "ZT16"))
zlm_model.ZT12 <- zlm(~ZT + cdr + (1|sample), 
                 data = sca,
                 method = 'glmer',
                 ebayes = FALSE,
                 parallel = TRUE)


# perform contrasts ----
# de_ZT4_vs_ZT0 <- summary(zlm_model.ZT0, doLRT = "ZTZT4")
# de_ZT12_vs_ZT0 <- summary(zlm_model.ZT0, doLRT = "ZTZT12")
# de_ZT16_vs_ZT0 <- summary(zlm_model.ZT0, doLRT = "ZTZT16")
de_ZT12_vs_ZT4 <- summary(zlm_model.ZT4, doLRT = "ZTZT12")
de_ZT16_vs_ZT4 <- summary(zlm_model.ZT4, doLRT = "ZTZT16")
de_ZT16_vs_ZT12 <- summary(zlm_model.ZT12, doLRT = "ZTZT16")

log2FC_threshold <- 0.585

# ZT4 vs ZT0
dt_ZT4_vs_ZT0 <- de_ZT4_vs_ZT0$datatable

results_ZT4_vs_ZT0 <- merge(
  dt_ZT4_vs_ZT0[contrast == "ZTZT4" & component == "H", .(primerid, `Pr(>Chisq)`)],
  dt_ZT4_vs_ZT0[contrast == "ZTZT4" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
  by = "primerid"
)

results_ZT4_vs_ZT0[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
results_ZT4_vs_ZT0[, log2FC := log2(exp(coef))]
results_ZT4_vs_ZT0 <- results_ZT4_vs_ZT0 |> 
  rename(gene = primerid) |>
  relocate(log2FC, fdr, .after = gene)

results_ZT4_vs_ZT0_signif <- results_ZT4_vs_ZT0 |> 
  filter(fdr < 0.05, abs(log2FC) > log2FC_threshold) |> 
  arrange(desc(log2FC)) |> 
  print()


# ZT12 vs ZT0
dt_ZT12_vs_ZT0 <- de_ZT12_vs_ZT0$datatable

results_ZT12_vs_ZT0 <- merge(
  dt_ZT12_vs_ZT0[contrast == "ZTZT12" & component == "H", .(primerid, `Pr(>Chisq)`)],
  dt_ZT12_vs_ZT0[contrast == "ZTZT12" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
  by = "primerid"
)

results_ZT12_vs_ZT0[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
results_ZT12_vs_ZT0[, log2FC := log2(exp(coef))]
results_ZT12_vs_ZT0 <- results_ZT12_vs_ZT0 |> 
  rename(gene = primerid) |>
  relocate(log2FC, fdr, .after = gene)

results_ZT12_vs_ZT0_signif <- results_ZT12_vs_ZT0 |> 
  filter(fdr < 0.05, abs(log2FC) > log2FC_threshold) |> 
  arrange(desc(log2FC)) |> 
  print()


# ZT16 vs ZT0
dt_ZT16_vs_ZT0 <- de_ZT16_vs_ZT0$datatable

results_ZT16_vs_ZT0 <- merge(
  dt_ZT16_vs_ZT0[contrast == "ZTZT16" & component == "H", .(primerid, `Pr(>Chisq)`)],
  dt_ZT16_vs_ZT0[contrast == "ZTZT16" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
  by = "primerid"
)

results_ZT16_vs_ZT0[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
results_ZT16_vs_ZT0[, log2FC := log2(exp(coef))]
results_ZT16_vs_ZT0 <- results_ZT16_vs_ZT0 |> 
  rename(gene = primerid) |>
  relocate(log2FC, fdr, .after = gene)

results_ZT16_vs_ZT0_signif <- results_ZT16_vs_ZT0 |> 
  filter(fdr < 0.05, abs(log2FC) > log2FC_threshold) |> 
  arrange(desc(log2FC)) |> 
  print()


# ZT12 vs ZT4
dt_ZT12_vs_ZT4 <- de_ZT12_vs_ZT4$datatable

results_ZT12_vs_ZT4 <- merge(
  dt_ZT12_vs_ZT4[contrast == "ZTZT12" & component == "H", .(primerid, `Pr(>Chisq)`)],
  dt_ZT12_vs_ZT4[contrast == "ZTZT12" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
  by = "primerid"
)

results_ZT12_vs_ZT4[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
results_ZT12_vs_ZT4[, log2FC := log2(exp(coef))]
results_ZT12_vs_ZT4 <- results_ZT12_vs_ZT4 |> 
  rename(gene = primerid) |>
  relocate(log2FC, fdr, .after = gene)

results_ZT12_vs_ZT4_signif <- results_ZT12_vs_ZT4 |> 
  filter(fdr < 0.05, abs(log2FC) > log2FC_threshold) |> 
  arrange(desc(log2FC)) |> 
  print()


# ZT16 vs ZT4
dt_ZT16_vs_ZT4 <- de_ZT16_vs_ZT4$datatable

results_ZT16_vs_ZT4 <- merge(
  dt_ZT16_vs_ZT4[contrast == "ZTZT16" & component == "H", .(primerid, `Pr(>Chisq)`)],
  dt_ZT16_vs_ZT4[contrast == "ZTZT16" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
  by = "primerid"
)

results_ZT16_vs_ZT4[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
results_ZT16_vs_ZT4[, log2FC := log2(exp(coef))]
results_ZT16_vs_ZT4 <- results_ZT16_vs_ZT4 |> 
  rename(gene = primerid) |>
  relocate(log2FC, fdr, .after = gene)

results_ZT16_vs_ZT4_signif <- results_ZT16_vs_ZT4 |> 
  filter(fdr < 0.05, abs(log2FC) > log2FC_threshold) |> 
  arrange(desc(log2FC)) |> 
  print()


# ZT16 vs ZT12
dt_ZT16_vs_ZT12 <- de_ZT16_vs_ZT12$datatable

results_ZT16_vs_ZT12 <- merge(
  dt_ZT16_vs_ZT12[contrast == "ZTZT16" & component == "H", .(primerid, `Pr(>Chisq)`)],
  dt_ZT16_vs_ZT12[contrast == "ZTZT16" & component == "logFC", .(primerid, coef, ci.hi, ci.lo)],
  by = "primerid"
)

results_ZT16_vs_ZT12[, fdr := p.adjust(`Pr(>Chisq)`, 'fdr')]
results_ZT16_vs_ZT12[, log2FC := log2(exp(coef))]
results_ZT16_vs_ZT12 <- results_ZT16_vs_ZT12 |> 
  rename(gene = primerid) |>
  relocate(log2FC, fdr, .after = gene)

results_ZT16_vs_ZT12_signif <- results_ZT16_vs_ZT12 |> 
  filter(fdr < 0.05, abs(log2FC) > log2FC_threshold) |> 
  arrange(desc(log2FC)) |> 
  print()


# save results in files
write_csv(results_ZT4_vs_ZT0, "04-analysis/DEGs/Dec2024_ZT/SE_ZT4_vs_ZT0.csv")
write_csv(results_ZT12_vs_ZT0, "04-analysis/DEGs/Dec2024_ZT/SE_ZT12_vs_ZT0.csv")
write_csv(results_ZT16_vs_ZT0, "04-analysis/DEGs/Dec2024_ZT/SE_ZT16_vs_ZT0.csv")
write_csv(results_ZT12_vs_ZT4, "04-analysis/DEGs/Dec2024_ZT/SE_ZT12_vs_ZT4.csv")
write_csv(results_ZT16_vs_ZT4, "04-analysis/DEGs/Dec2024_ZT/SE_ZT16_vs_ZT4.csv")
write_csv(results_ZT16_vs_ZT12, "04-analysis/DEGs/Dec2024_ZT/SE_ZT16_vs_ZT12.csv")

# heatmap ----
# assemble all genes and their contrast of origin
results_signif_all <- bind_rows(
  results_ZT4_vs_ZT0_signif |> mutate(contrast = "ZT4_vs_ZT0"),
  results_ZT12_vs_ZT0_signif |> mutate(contrast = "ZT12_vs_ZT0"),
  results_ZT16_vs_ZT0_signif |> mutate(contrast = "ZT16_vs_ZT0"),
  results_ZT12_vs_ZT4_signif |> mutate(contrast = "ZT12_vs_ZT4"),
  results_ZT16_vs_ZT4_signif |> mutate(contrast = "ZT16_vs_ZT4"),
  results_ZT16_vs_ZT12_signif |> mutate(contrast = "ZT16_vs_ZT12")
) |> 
  relocate(contrast, .after = gene) |> 
  print()
write_csv(results_signif_all, "04-analysis/DEGs/Dec2024_ZT/SE_all_signif_genes.csv")

# remove duplicates
zt_genes <- results_signif_all |> 
  distinct(gene) |> 
  pull()

# make a "collection time" for samples
nuc.subset <- nuclei |> 
  subset(subclass_name == '016 CA1-ProS Glut') |>
  subset(activity_condition == 'SE' | activity_condition == 'EE6h')

# for cases where activity_condition == EE6h, set ZT_collection to ZT + 6. Otherwise, use ZT. Please use case_when
nuc.subset@meta.data <- nuc.subset@meta.data |> 
  mutate(ZT.collection = case_when(
    (ZT == 'ZT0' & activity_condition == 'EE6h') ~ 'ZT6',
    (ZT == 'ZT4' & activity_condition == 'EE6h') ~ 'ZT10',
    (ZT == 'ZT12' & activity_condition == 'EE6h') ~ 'ZT18',
    (ZT == 'ZT16' & activity_condition == 'EE6h') ~ 'ZT22',
    TRUE ~ ZT
  )) |> 
  mutate(ZT.collection = factor(ZT.collection, 
    levels = c('ZT0', 'ZT4', 'ZT6', 'ZT10', 'ZT12', 'ZT16', 'ZT18', 'ZT22')))
  
# get average expression for each ZT gene
# mat_zt <- AverageExpression(nuc.subset, group.by = 'ZT.collection', features = zt_genes)$SCT |> 
mat_zt <- AverageExpression(celltype, group.by = 'ZT', features = zt_genes)$SCT |> 
  t() |> 
  # z-score the genes
  scale()


k = 4
col_fun <- circlize::colorRamp2(c(-2,0,2), hcl_palette = "blue-red2")
hm <- Heatmap(t(mat_zt),
              # col = col_fun,
              cluster_columns = FALSE,
              row_split = k,
              cluster_row_slices = F
              # row_names_gp = gpar(fontsize = 5)
              )
hm = draw(hm)

InteractiveComplexHeatmap::htShiny(hm)


if (SAVE_PLOTS) {
  png("05-results/figureZT/raw_R_plots/SE_ZT_DEGs.png", width = 5, height = 10, units = "in", res = 900)
  print(hm)
  dev.off()
}

# GO analysis ----
dend <- row_dend(hm)
row_clusters <- cutree(as.hclust(dend), k = k)

# run GO on each set
go_terms_all <- tibble()

terms.1 <- RunGOEnrichment(names(row_clusters[row_clusters == 1])) |> mutate(cluster = 1)
go_terms_all <- rbind(go_terms_all, terms.1)

terms.2 <- RunGOEnrichment(names(row_clusters[row_clusters == 2])) |> mutate(cluster = 2)
go_terms_all <- rbind(go_terms_all, terms.2)

terms.3 <- RunGOEnrichment(names(row_clusters[row_clusters == 3])) |> mutate(cluster = 3)
go_terms_all <- rbind(go_terms_all, terms.3)

terms.4 <- RunGOEnrichment(names(row_clusters[row_clusters == 4])) |> mutate(cluster = 4)
go_terms_all <- rbind(go_terms_all, terms.4)

terms.5 <- RunGOEnrichment(names(row_clusters[row_clusters == 5])) |> mutate(cluster = 5)
go_terms_all <- rbind(go_terms_all, terms.5)