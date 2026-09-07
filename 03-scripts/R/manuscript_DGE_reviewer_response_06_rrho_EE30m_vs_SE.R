# Reviewer response follow-up: rank-rank hypergeometric overlap (RRHO)
# comparing the DESeq2 pseudobulk gene ranking against two alternative
# methods -- limma-voom pseudobulk and a negative-binomial GLMM with an
# animal random effect -- for 016 CA1-ProS Glut, EE30m vs. SE.
#
# Styled after Figure 1A/1D of Piron et al. 2024, Life Science Alliance
# (RedRibbon paper, https://www.life-science-alliance.org/content/7/2/e202302203):
# each comparison gets a "real data" map, which should show a diagonal signal
# from down- to up-regulation when the two rankings agree (their Fig 1A), and
# a negative-control map with one list's gene order scrambled, which should
# show no signal (their Fig 1D).
#
# manuscript_DGE_reviewer_response_07_rrho_EE6h_vs_SE.R is the same analysis
# at the 6 h timepoint. Between them the two scripts produce the four panels
# the reviewer response needs, all on one gene set at matched N.
#
# WHICH GENES
#
# All four comparisons run on a single pinned gene set:
#   04-analysis/reviewer_response/rrho/016_CA1-ProS_Glut__analysis_genes_stratified.csv
# 2015 genes drawn equally from every DESeq2 baseMean decile, plus the
# 15-gene IEG panel. Fixing one set across both timepoints and both methods is
# what makes the four maps comparable to each other: -log10(p) on an RRHO map
# scales with the number of genes, so maps at different N cannot be read
# against one another even when the rankings agree identically well.
#
# The NB-GLMM has to be fit per gene (glmmTMB, ~2.8 sec/gene, i.e. ~15 hours
# serial genome-wide for one contrast), so a subset is unavoidable and the
# question is which subset. An earlier version of this script took the top
# 2000 genes by mean per-cell SCT expression, which was a poor choice for a
# method comparison: it covered only the top ~11% of the abundance range,
# contained none of the 15 IEGs (the per-cell mean rewards broadly-detected
# genes over sparsely-but-strongly induced ones), and landed on the easiest
# genes -- DESeq2-vs-voom concordance here runs 0.997 in the top abundance
# decile but 0.955 in the bottom, so a top-N set overstates agreement (0.997
# vs 0.986 genome-wide). The stratified set pools to 0.985, essentially the
# genome-wide value. See
# 06-reports/exploratory_notebooks/investigate_RRHO_map_saturation.qmd and
# SelectStratifiedGenes()/ThinStratifiedGenes() in seq_functions.R.
#
# RRHO method: genes are ranked within each list by a signed statistic
# (sign(log2FC) * -log10(p), RRHORankStat()), from most down-regulated to
# most up-regulated. At a grid of rank-cutoff pairs (i, j), the number of
# genes in common between the top-i of list 1 and top-j of list 2 is tested
# against a hypergeometric null; the resulting -log10(p) grid is the RRHO
# map. This is the test Plaisier et al. 2010 introduced and Piron et al. 2024
# describes in their Methods. computeRRHO() takes that tail on the log scale,
# which is load-bearing at these list lengths.

library(tidyverse)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/reviewer_response"
deseq_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
voom_dir <- file.path(out_dir, "voom")
rrho_dir <- file.path(out_dir, "rrho")
dir.create(rrho_dir, showWarnings = FALSE, recursive = TRUE)

subclass <- "016 CA1-ProS Glut"
subclass_fname <- ShrinkSubclassName(subclass)
group1 <- "EE30m"
group2 <- "SE"

stem <- glue("{subclass_fname}__{group1}_vs_{group2}")
genes_path <- file.path(rrho_dir, glue("{subclass_fname}__analysis_genes_stratified.csv"))


# ============================================================================ #
# gene set ----
# ============================================================================ #
if (!file.exists(genes_path)) {
  stop(glue("{genes_path} not found. It is the abundance-stratified set shared ",
            "by all four RRHO comparisons; see the header for how it was built."))
}
gene_set <- read_csv(genes_path, show_col_types = FALSE)
analysis_genes <- gene_set$gene
message(glue("Gene set: {length(analysis_genes)} genes across ",
             "{n_distinct(gene_set$stratum)} abundance strata ",
             "({sum(gene_set$forced)} force-included IEGs)"))


# ============================================================================ #
# rankings ----
# ============================================================================ #
deseq <- read_csv(file.path(deseq_dir, glue("{stem}.csv")), show_col_types = FALSE) |>
  filter(gene %in% analysis_genes) |>
  drop_na(log2FoldChange.shrink, pvalue) |>
  transmute(gene,
            log2fc = log2FoldChange.shrink,
            signed_stat = RRHORankStat(log2FoldChange.shrink, pvalue))

voom <- read_csv(file.path(voom_dir, glue("{stem}.csv")), show_col_types = FALSE) |>
  filter(gene %in% analysis_genes) |>
  drop_na(log2FoldChange, pvalue) |>
  transmute(gene,
            log2fc = log2FoldChange,
            signed_stat = RRHORankStat(log2FoldChange, pvalue))

glmm <- FitSubclassGLMM(
  analysis_genes, subclass_fname, group1, group2,
  cache_path = file.path(rrho_dir, glue("{stem}__NB-GLMM_stratified.csv"))
) |>
  drop_na(log2FoldChange, pvalue) |>
  transmute(gene,
            log2fc = log2FoldChange,
            signed_stat = RRHORankStat(log2FoldChange, pvalue))

message(glue("Ranked genes -- DESeq2: {nrow(deseq)}, limma-voom: {nrow(voom)}, NB-GLMM: {nrow(glmm)}"))


# ============================================================================ #
# the two comparisons ----
# ============================================================================ #
# RRHOCompare() writes the map, the scrambled negative control, and both
# quadrant tables, and returns the plots plus a concordance row.
res_voom <- RRHOCompare(deseq, voom, "limma-voom",
                        out_stem = glue("{stem}__voom"), rrho_dir = rrho_dir)
res_glmm <- RRHOCompare(deseq, glmm, "NB-GLMM",
                        out_stem = glue("{stem}__glmm"), rrho_dir = rrho_dir)

concordance <- bind_rows(res_voom$concordance, res_glmm$concordance) |>
  mutate(contrast = glue("{group1} vs {group2}"), .before = 1)
write_csv(concordance, file.path(rrho_dir, glue("{stem}__rank_concordance.csv")))
print(concordance)

quadrants <- bind_rows(
  res_voom$quadrants |> mutate(comparison = "DESeq2 ~ limma-voom", control = FALSE, .before = 1),
  res_glmm$quadrants |> mutate(comparison = "DESeq2 ~ NB-GLMM", control = FALSE, .before = 1),
  res_voom$quadrants_scrambled |> mutate(comparison = "DESeq2 ~ limma-voom", control = TRUE, .before = 1),
  res_glmm$quadrants_scrambled |> mutate(comparison = "DESeq2 ~ NB-GLMM", control = TRUE, .before = 1)
) |>
  mutate(contrast = glue("{group1} vs {group2}"), .before = 1)
write_csv(quadrants, file.path(rrho_dir, glue("{stem}__RRHO_quadrants_all.csv")))


# ============================================================================ #
# figures ----
# ============================================================================ #
# One shared colour scale across the two real maps, so the panels can be read
# against each other rather than each against its own maximum.
shared_cap <- max(quantile(res_voom$map$neglog10p, 0.99),
                  quantile(res_glmm$map$neglog10p, 0.99))

p_main <-
  plotRRHO(res_voom$map, "DESeq2 vs. limma-voom",
          y_label = "limma-voom rank (down \u2192 up)", cap = shared_cap) +
    # theme(axis.title.y = element_text(family = "Arial")) +
plotRRHO(res_glmm$map, "DESeq2 vs. NB-GLMM",
           y_label = "NB-GLMM rank (down \u2192 up)", cap = shared_cap) +
  # theme(axis.title.y = element_text(family = "Arial")) +
  patchwork::plot_annotation(
    title = glue("{subclass}: {group1} vs {group2}")
  )
print(p_main)

p_controls <-
  res_voom$plot_scrambled + res_glmm$plot_scrambled +
  patchwork::plot_annotation(
    title = glue("RRHO negative controls ({subclass}, {group1} vs {group2})")
  )
print(p_controls)

if (exists("SAVE_PLOTS") && SAVE_PLOTS) {
  for (ext in c("png", "svg")) {
    ggsave(file.path(rrho_dir, glue("RRHO_{group1}_vs_{group2}.{ext}")), p_main, width = 11, height = 5)
    ggsave(file.path(rrho_dir, glue("RRHO_{group1}_vs_{group2}_negative_controls.{ext}")), p_controls, width = 11, height = 5)
  }
}

