# Reviewer response follow-up: rank-rank hypergeometric overlap (RRHO)
# comparing DESeq2 pseudobulk vs. limma-voom pseudobulk gene rankings for one
# subclass x contrast, styled after Figure 1A/1D of Piron et al. 2024, Life
# Science Alliance (RedRibbon paper,
# https://www.life-science-alliance.org/content/7/2/e202302203):
#   - a "real data" map showing a diagonal signal from down- to
#     up-regulation when the two rankings agree (like their Fig 1A)
#   - a negative-control map with one list's gene order scrambled, which
#     should show no signal (like their Fig 1D)
#
# Subclass/contrast: 016 CA1-ProS Glut, EE30m vs. SE. Unlike
# manuscript_DGE_reviewer_response_06_rrho_deseq2_vs_glmm.R (which had to fit
# a genome-wide NB-GLMM from scratch), both DESeq2 and limma-voom results for
# this subclass x contrast are already genome-wide on disk, so this script
# needs no Seurat object and no model fitting.
#
# RRHO method: genes are ranked within each list by a signed statistic
# (sign(log2FC) * -log10(p)), from most down-regulated to most
# up-regulated. At a grid of rank-cutoff pairs (i, j), the number of genes
# in common between the top-i of list 1 and top-j of list 2 is tested
# against a hypergeometric null; the resulting -log10(p) grid is the RRHO
# map. computeRRHO()/plotRRHO() are shared helpers in seq_functions.R.

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


# ============================================================================ #
# signed rankings, DESeq2 vs. limma-voom, shared gene set ----
# ============================================================================ #
deseq <- read_csv(file.path(deseq_dir, glue("{subclass_fname}__{group1}_vs_{group2}.csv")),
                  show_col_types = FALSE) |>
  drop_na(log2FoldChange.shrink, pvalue) |>
  transmute(gene, signed_stat = sign(log2FoldChange.shrink) * -log10(pmax(pvalue, 1e-300)))

voom <- read_csv(file.path(voom_dir, glue("{subclass_fname}__{group1}_vs_{group2}.csv")),
                 show_col_types = FALSE) |>
  drop_na(log2FoldChange, pvalue) |>
  transmute(gene, signed_stat = sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300)))

x_label <- "DESeq2 rank (down → up)"
y_label <- "limma-voom rank (down → up)"

# --- real-data RRHO map, full shared gene set ----
df_rrho_real <- computeRRHO(deseq, voom)
write_csv(df_rrho_real, file.path(rrho_dir, glue("{subclass_fname}__{group1}_vs_{group2}__RRHO_voom_real.csv")))

p_real <- plotRRHO(df_rrho_real,
                   glue("RRHO: DESeq2 vs. limma-voom ({subclass}, {group1} vs {group2})"),
                   x_label = x_label, y_label = y_label)
ggsave(file.path(rrho_dir, "RRHO_voom_real.pdf"), p_real, width = 6, height = 5.5)

# --- negative control: scramble the voom gene order ----
set.seed(17)
voom_scrambled <- voom |>
  mutate(signed_stat = sample(signed_stat))

df_rrho_scrambled <- computeRRHO(deseq, voom_scrambled)
write_csv(df_rrho_scrambled, file.path(rrho_dir, glue("{subclass_fname}__{group1}_vs_{group2}__RRHO_voom_scrambled.csv")))

p_scrambled <- plotRRHO(df_rrho_scrambled,
                        glue("RRHO negative control: DESeq2 vs. scrambled limma-voom ({subclass}, {group1} vs {group2})"),
                        x_label = x_label, y_label = y_label)
ggsave(file.path(rrho_dir, "RRHO_voom_scrambled_negative_control.pdf"), p_scrambled, width = 6, height = 5.5)


# ============================================================================ #
# secondary map: same top-2000-gene set used by the DESeq2-vs-GLMM RRHO ----
# (directly comparable to manuscript_DGE_reviewer_response_06's map)
# ============================================================================ #
top2000_path <- file.path(rrho_dir, glue("{subclass_fname}__{group1}_vs_{group2}__NB-GLMM_top2000genes.csv"))
if (file.exists(top2000_path)) {
  top_genes <- read_csv(top2000_path, show_col_types = FALSE) |> pull(gene) |> unique()

  deseq_top <- deseq |> filter(gene %in% top_genes)
  voom_top  <- voom  |> filter(gene %in% top_genes)

  df_rrho_top2000 <- computeRRHO(deseq_top, voom_top)
  write_csv(df_rrho_top2000, file.path(rrho_dir, glue("{subclass_fname}__{group1}_vs_{group2}__RRHO_voom_top2000.csv")))

  p_top2000 <- plotRRHO(df_rrho_top2000,
                        glue("RRHO (top 2000 genes): DESeq2 vs. limma-voom ({subclass}, {group1} vs {group2})"),
                        x_label = x_label, y_label = y_label)
  ggsave(file.path(rrho_dir, "RRHO_voom_top2000.pdf"), p_top2000, width = 6, height = 5.5)
} else {
  message(glue("Skipping top-2000-gene RRHO map: {top2000_path} not found (run manuscript_DGE_reviewer_response_06_rrho_deseq2_vs_glmm.R first if you want this comparison)."))
}

print(glue("RRHO (DESeq2 vs. limma-voom) analysis complete. Outputs in {rrho_dir}/"))
