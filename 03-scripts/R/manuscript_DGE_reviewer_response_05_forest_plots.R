# Reviewer response (d): forest plots of effect sizes with 95% CIs, comparing
# methods, to improve interpretability beyond binary FDR thresholds.
# Run these first: 00_cache_subclasses.R, 02_pseudobulk_voom.R,
# 03_glmm_gene_expression.R
#
# For the IEG panel, per subclass x contrast: log2FC point estimates with 95%
# CIs from DESeq2 (apeglm-shrunken, CI = estimate +/- 1.96*lfcSE), limma-voom,
# and the NB-GLMM with animal as a random effect. A second plot set adds a
# 4th series: DESeq2's raw (unshrunk) Wald estimate, i.e. before apeglm --
# useful for seeing how much the shrinkage step itself moved the estimate,
# separate from cross-method disagreement.
#
# Some NB-GLMM fits (small subclasses / low cell counts per condition) don't
# converge to a stable estimate -- their CI can blow up to a width many
# orders of magnitude larger than the point estimate itself (e.g. a log2FC of
# 46 with a CI of +/-20 million), which flattens every other point on the same
# plot. Those are dropped before plotting: any NB-GLMM row with a non-finite
# CI, or with CI width >= 10x |log2FoldChange|, is excluded. That threshold
# comfortably clears normal estimation noise (median CI-width:effect ratio
# among GLMM fits is ~1.6, 3rd quartile ~4.4) while catching the clear
# separation-driven blowups (ratios up to ~9x10^5 in this data).

library(tidyverse)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/reviewer_response"
deseq_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
dir.create(file.path(out_dir, "forest_plots"), showWarnings = FALSE, recursive = TRUE)

subclass_list <- LoadSubclassesToUse()
log2FC_threshold <- 0.585
ieg_panel <- LoadGeneList("IEG")

glmm_ci_ratio_max <- 10  # exclude fits where CI width >= this x |log2FoldChange|

df_glmm_all <- read_csv(file.path(out_dir, "glmm", "IEG_NB-GLMM_animal-random-effect.csv"),
                        show_col_types = FALSE) |>
  filter(
    is.finite(log2FC_CI_low), is.finite(log2FC_CI_high),
    (log2FC_CI_high - log2FC_CI_low) < glmm_ci_ratio_max * pmax(abs(log2FoldChange), 1e-6)
  )


# ============================================================================ #
# IEG log2FC forest plots: DESeq2 vs. limma-voom vs. NB-GLMM ----
# ============================================================================ #
collectForestData <- function(subclass) {
  subclass_fname <- ShrinkSubclassName(subclass)
  contrast_files <- list.files(deseq_dir, pattern = glue("^{subclass_fname}__"))

  map_dfr(contrast_files, function(f) {
    ct <- str_remove(f, glue("^{subclass_fname}__")) |> str_remove(".csv$")

    deseq_raw_shrink <- read_csv(file.path(deseq_dir, f), show_col_types = FALSE) |>
      filter(gene %in% ieg_panel)

    deseq <- deseq_raw_shrink |>
      transmute(gene,
                method = "DESeq2 (apeglm)",
                log2FoldChange = log2FoldChange.shrink,
                CI_low  = log2FoldChange.shrink - 1.96 * lfcSE.shrink,
                CI_high = log2FoldChange.shrink + 1.96 * lfcSE.shrink)

    deseq_raw <- deseq_raw_shrink |>
      transmute(gene,
                method = "DESeq2 (raw, unshrunk)",
                log2FoldChange = log2FoldChange.raw,
                CI_low  = log2FoldChange.raw - 1.96 * lfcSE.raw,
                CI_high = log2FoldChange.raw + 1.96 * lfcSE.raw)

    voom_path <- file.path(out_dir, "voom", f)
    voom <- if (file.exists(voom_path)) {
      read_csv(voom_path, show_col_types = FALSE) |>
        filter(gene %in% ieg_panel) |>
        transmute(gene, method = "limma-voom",
                  log2FoldChange, CI_low = log2FC_CI_low, CI_high = log2FC_CI_high)
    } else tibble()

    glmm <- df_glmm_all |>
      filter(subclass == !!subclass, contrast == ct) |>
      transmute(gene, method = "NB-GLMM (animal RE)",
                log2FoldChange, CI_low = log2FC_CI_low, CI_high = log2FC_CI_high)

    bind_rows(deseq, deseq_raw, voom, glmm) |>
      mutate(subclass = subclass, contrast = ct)
  })
}

df_forest_4method <- map_dfr(subclass_list, collectForestData)
write_csv(df_forest_4method, file.path(out_dir, "0_forest_plot_data.csv"))

df_forest_3method <- df_forest_4method |> filter(method != "DESeq2 (raw, unshrunk)")

method_colors_3 <- c("DESeq2 (apeglm)" = "#1E88E5",
                     "limma-voom" = "#D81B60",
                     "NB-GLMM (animal RE)" = "#FFC107")

method_colors_4 <- c(method_colors_3, "DESeq2 (raw, unshrunk)" = "#098154")

plotForestSet <- function(df_forest, method_colors, pdf_path) {
  pdf(pdf_path, width = 8, height = 6)
  for (sc in unique(df_forest$subclass)) {
    df_sc <- df_forest |> filter(subclass == sc)
    for (ct in unique(df_sc$contrast)) {
      df_plot <- df_sc |>
        filter(contrast == ct) |>
        mutate(gene = fct_rev(factor(gene, levels = sort(ieg_panel))))
      if (nrow(df_plot) == 0) next

      p <- ggplot(df_plot,
                  aes(x = log2FoldChange, y = gene, color = method)) +
        geom_vline(xintercept = 0, linetype = "solid", color = "gray60") +
        geom_vline(xintercept = c(-log2FC_threshold, log2FC_threshold),
                   linetype = "dashed", color = "gray80") +
        geom_pointrange(aes(xmin = CI_low, xmax = CI_high),
                        position = position_dodge(width = 0.6),
                        size = 0.35) +
        scale_color_manual(values = method_colors) +
        labs(title = glue("{sc}: {str_replace(ct, '_vs_', ' vs. ')}"),
             x = "log2 fold-change (95% CI)", y = NULL, color = NULL) +
        theme_minimal(base_family = "Helvetica") +
        theme(legend.position = "bottom")
      print(p)
    }
  }
  dev.off()
}

plotForestSet(df_forest_3method, method_colors_3,
             file.path(out_dir, "forest_plots", "IEG_forest_plots.pdf"))
plotForestSet(df_forest_4method, method_colors_4,
             file.path(out_dir, "forest_plots", "IEG_forest_plots_with_raw_deseq2.pdf"))

print(glue("Forest plots complete: {out_dir}/forest_plots/"))
