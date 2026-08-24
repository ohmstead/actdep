# Reviewer response (a): per-subclass effective sample sizes & replicate numbers.
# Run manuscript_DGE_reviewer_response_00_cache_subclasses.R first.
#
# For each subclass x condition we report:
#   n_animals   : biological replicates contributing >= 1 nucleus
#   n_cells     : nuclei
#   cells/animal: min / median / max
#   ICC         : intra-class correlation (between-animal variance share),
#                 median across the expressed IEG panel, estimated with
#                 lmer(expr ~ activity_condition + (1|sample))
#   n_eff       : n_cells / (1 + (m_bar - 1) * ICC), the design-effect-adjusted
#                 effective number of independent nuclei (Kish 1965)
#
# NOTE: each sample is one animal (one hemisphere per animal was sequenced),
# so `sample` is the animal-level unit and (1|sample) in the ICC model is a
# genuine animal-level random effect.

library(Seurat)
library(tidyverse)
library(lme4)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/reviewer_response"
cache_dir <- file.path(out_dir, "subclass_cache")

subclass_list <- LoadSubclassesToUse()
ieg_panel <- LoadGeneList("IEG")

estimateSubclassICC <- function(seurat_subclass, genes) {
  expr_mat <- GetAssayData(seurat_subclass, assay = "SCT", layer = "data")
  genes <- intersect(genes, rownames(expr_mat))
  genes <- genes[rowMeans(expr_mat[genes, , drop = FALSE] > 0) >= 0.01]
  if (length(genes) == 0) return(NA_real_)

  cell_df <- seurat_subclass@meta.data |>
    as_tibble(rownames = "cell") |>
    select(cell, sample, activity_condition)

  iccs <- map_dbl(genes, function(g) {
    cell_df$expr <- expr_mat[g, cell_df$cell]
    fit <- tryCatch(
      lmer(expr ~ activity_condition + (1 | sample), data = cell_df,
           control = lmerControl(calc.derivs = FALSE)),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    vc <- as.data.frame(VarCorr(fit))
    var_animal <- vc$vcov[vc$grp == "sample"]
    var_resid  <- vc$vcov[vc$grp == "Residual"]
    var_animal / (var_animal + var_resid)
  })

  median(iccs, na.rm = TRUE)
}

computeESS <- function(nuclei_subclass, subclass) {
  icc <- estimateSubclassICC(nuclei_subclass, ieg_panel)

  nuclei_subclass@meta.data |>
    as_tibble() |>
    group_by(activity_condition, sample) |>
    summarize(n_cells = n(), .groups = "drop") |>
    group_by(activity_condition) |>
    summarize(
      n_animals = n_distinct(sample),
      cells_per_animal_min = min(n_cells),
      cells_per_animal_median = median(n_cells),
      cells_per_animal_max = max(n_cells),
      n_cells = sum(n_cells),
      .groups = "drop"
    ) |>
    mutate(
      subclass = subclass,
      ICC_IEG_median = icc,
      mean_cells_per_animal = n_cells / n_animals,
      design_effect = 1 + (mean_cells_per_animal - 1) * icc,
      n_eff_cells = n_cells / design_effect
    ) |>
    relocate(subclass)
}

df_ess <- map_dfr(subclass_list, function(subclass) {
  message(glue("ESS: {subclass}"))
  subclass_fname <- ShrinkSubclassName(subclass)
  nuclei_subclass <- readRDS(file.path(cache_dir, glue("{subclass_fname}.Rds")))
  res <- computeESS(nuclei_subclass, subclass)
  rm(nuclei_subclass)
  gc()
  res
})

write_csv(df_ess, file.path(out_dir, "0_per_subclass_replicates_and_ESS.csv"))
print(df_ess, n = Inf)
print(glue("Sample-size report complete: {out_dir}/0_per_subclass_replicates_and_ESS.csv"))
