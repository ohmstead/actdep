# Reviewer response (d), CA1 panel: forest plot of IEG effect sizes with 95%
# CIs across all four DGE methods, in a single 2x2 contrast grid.
# Run these first: 00_cache_subclasses.R, 02_pseudobulk_voom.R,
# 03_glmm_gene_expression.R
#
# Companion to 05_forest_plots.R, which loops over every subclass x contrast
# and emits one page per combination. This script instead pins the subclass to
# 016 CA1-ProS Glut and lays the four stimulus x timepoint contrasts out as a
# 2x2 grid on one page (columns = stimulus EE/KA, rows = timepoint 30m/6h), so
# the reader can compare method agreement across contrasts at a glance.
#
# Only the 15 IEGs are plotted -- circadian genes are excluded here (the
# original script's IEG + circadian union is kept there). All four estimates
# are shown per gene: DESeq2's raw (unshrunk) Wald estimate, the apeglm-shrunken
# estimate, limma-voom, and the NB-GLMM with animal as a random effect. CIs are
# estimate +/- 1.96*SE for the DESeq2 series; voom and GLMM carry their own.
#
# NB-GLMM fits that failed to converge are dropped before plotting, but the
# criterion differs from 05_forest_plots.R. That script drops a fit when the CI
# width is >= 10x |log2FoldChange|, which works when scanning every subclass
# (including small ones whose fits diverge to CIs of +/- millions). Applied to
# just the 15 IEGs in CA1 it misfires: a gene with a genuinely null effect has a
# tiny |log2FoldChange| in the denominator, so a perfectly healthy CI of width
# ~1 yields a huge ratio and gets cut. That silently removes the null GLMM
# results and makes the method look more concordant than it is. Here the cut is
# an absolute CI width instead, which still catches real divergence without
# penalising near-zero effects. Across all CA1 IEG fits the widest CI is ~5.2,
# so nothing is currently dropped -- the filter is purely defensive.

library(tidyverse)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/reviewer_response"
deseq_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
plot_dir <- file.path(out_dir, "forest_plots")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

subclass <- "016 CA1-ProS Glut"
subclass_label <- "CA1-ProS Glut"
log2FC_threshold <- 0.585
ieg_genes <- LoadGeneList("IEG")
gene_order <- sort(ieg_genes)

# 2x2 grid, filled by row: top row = the 30 min contrasts, bottom = 6 h, with
# EE in the left column and KA in the right. Each panel is titled with its own
# full contrast name; stimulus/timepoint are kept as columns because the x-axis
# range is shared down each stimulus column (see df_axis_range below).
contrast_grid <- tribble(
  ~contrast,        ~contrast_label,      ~stimulus,     ~timepoint,
  "EE30m_vs_SE",    "EE30m vs. SE",       "EE",          "30 min",
  "KA30m_vs_SE",    "KA30m vs. SE",       "KA",          "30 min",
  "EE6h_vs_SE",     "EE6h vs. SE",        "EE",          "6 h",
  "KA6h_vs_SE",     "KA6h vs. SE",        "KA",          "6 h"
)
contrast_labels <- contrast_grid$contrast_label

method_levels <- c("DESeq2 (raw, unshrunk)", "DESeq2 (apeglm)",
                   "limma-voom", "NB-GLMM (animal RE)")
method_colors <- c("DESeq2 (raw, unshrunk)" = "gray30",
                   "DESeq2 (apeglm)"        = "black",
                   "limma-voom"             = "gray55",
                   "NB-GLMM (animal RE)"    = "gray75")

glmm_ci_width_max <- 20  # exclude fits whose 95% CI is wider than this (log2 units)

df_glmm_raw <- read_csv(file.path(out_dir, "glmm", "IEG_NB-GLMM_animal-random-effect.csv"),
                        show_col_types = FALSE) |>
  filter(subclass == !!subclass, gene %in% ieg_genes)

df_glmm_all <- df_glmm_raw |>
  filter(
    is.finite(log2FC_CI_low), is.finite(log2FC_CI_high),
    (log2FC_CI_high - log2FC_CI_low) < glmm_ci_width_max
  )

# Say out loud which GLMM fits the filter removed -- a silently dropped series
# leaves a gap in the forest plot that reads as a missing estimate.
n_dropped <- nrow(df_glmm_raw) - nrow(df_glmm_all)
if (n_dropped > 0) {
  df_glmm_raw |>
    anti_join(df_glmm_all, by = c("gene", "contrast")) |>
    transmute(gene, contrast, log2FoldChange,
              CI_width = log2FC_CI_high - log2FC_CI_low) |>
    print()
}
print(glue("NB-GLMM fits dropped as non-convergent: {n_dropped} of {nrow(df_glmm_raw)}"))


# ============================================================================ #
# Collect the four methods' estimates for each contrast ----
# ============================================================================ #
collectForestData <- function(ct) {
  subclass_fname <- ShrinkSubclassName(subclass)
  f <- glue("{subclass_fname}__{ct}.csv")

  deseq_raw_shrink <- read_csv(file.path(deseq_dir, f), show_col_types = FALSE) |>
    filter(gene %in% ieg_genes)

  deseq_raw <- deseq_raw_shrink |>
    transmute(gene,
              method = "DESeq2 (raw, unshrunk)",
              log2FoldChange = log2FoldChange.raw,
              CI_low  = log2FoldChange.raw - 1.96 * lfcSE.raw,
              CI_high = log2FoldChange.raw + 1.96 * lfcSE.raw)

  deseq <- deseq_raw_shrink |>
    transmute(gene,
              method = "DESeq2 (apeglm)",
              log2FoldChange = log2FoldChange.shrink,
              CI_low  = log2FoldChange.shrink - 1.96 * lfcSE.shrink,
              CI_high = log2FoldChange.shrink + 1.96 * lfcSE.shrink)

  voom <- read_csv(file.path(out_dir, "voom", f), show_col_types = FALSE) |>
    filter(gene %in% ieg_genes) |>
    transmute(gene, method = "limma-voom",
              log2FoldChange, CI_low = log2FC_CI_low, CI_high = log2FC_CI_high)

  glmm <- df_glmm_all |>
    filter(contrast == ct) |>
    transmute(gene, method = "NB-GLMM (animal RE)",
              log2FoldChange, CI_low = log2FC_CI_low, CI_high = log2FC_CI_high)

  bind_rows(deseq_raw, deseq, voom, glmm) |>
    mutate(subclass = subclass, contrast = ct)
}

df_forest <- map_dfr(contrast_grid$contrast, collectForestData) |>
  left_join(contrast_grid, by = "contrast") |>
  mutate(gene = fct_rev(factor(gene, levels = gene_order)),
         method = factor(method, levels = method_levels),
         contrast_label = factor(contrast_label, levels = contrast_labels),
         stimulus = factor(stimulus, levels = c("EE", "KA")),
         timepoint = factor(timepoint, levels = c("30 min", "6 h")))

write_csv(df_forest, file.path(plot_dir, "0_forest_plot_data_CA1.csv"))


# ============================================================================ #
# 2x2 forest plot ----
# ============================================================================ #
# facet_wrap gives each panel its own title but also its own x scale, which
# would break the 30 min vs 6 h comparison within a stimulus. Pad both panels of
# a stimulus with invisible points at that stimulus's overall min/max so the two
# panels in a column end up on a common range, as facet_grid gave us before.
df_axis_range <- df_forest |>
  group_by(stimulus) |>
  summarise(xmin = min(CI_low, na.rm = TRUE),
            xmax = max(CI_high, na.rm = TRUE), .groups = "drop") |>
  inner_join(contrast_grid, by = "stimulus", relationship = "one-to-many") |>
  pivot_longer(c(xmin, xmax), names_to = NULL, values_to = "log2FoldChange") |>
  mutate(contrast_label = factor(contrast_label, levels = contrast_labels),
         gene = df_forest$gene[1])

p <- ggplot(df_forest, aes(x = log2FoldChange, y = gene, color = method)) +
  geom_blank(data = df_axis_range, aes(x = log2FoldChange, y = gene),
             inherit.aes = FALSE) +
  geom_vline(xintercept = 0, linetype = "solid", color = "gray60") +
  geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold),
             linetype = "dashed", color = "gray80") +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
                  position = position_dodge(width = 0.7),
                  size = 0.3, fatten = 2) +
  scale_color_manual(values = method_colors, breaks = method_levels) +
  # Pin the gene axis explicitly. The geom_blank layer contains a single gene,
  # so letting ggplot train the discrete scale on layer order would put that one
  # gene first and append the remaining 14 after it.
  scale_y_discrete(limits = rev(gene_order)) +
  facet_wrap(~ contrast_label, nrow = 2, scales = "free_x", axes = "all_y") +
  labs(title = glue("{subclass_label}: IEG log2(FC) by method"),
       x = "log2 fold-change (95% CI)", y = NULL, color = NULL) +
  theme_minimal(base_family = "Helvetica") +
  theme(legend.position = "bottom",
        panel.spacing = unit(1, "lines"),
        strip.text = element_text(size = 12),
        axis.text.y = element_text(size = 8, face = "italic"))
print(p)

if (exists("SAVE_PLOTS") && SAVE_PLOTS) {
  ggsave(file.path(plot_dir, "IEG_forest_plots_CA1_2x2.svg"), p,
         width = 10, height = 9, bg = "white")
  ggsave(file.path(plot_dir, "IEG_forest_plots_CA1_2x2.png"), p,
         width = 10, height = 9, dpi = 300, bg = "white")
}
