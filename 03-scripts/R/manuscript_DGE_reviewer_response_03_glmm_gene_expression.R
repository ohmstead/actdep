# Reviewer response (c1): sample-level mixed-effects model for gene expression,
# animal as a random effect.
# Run manuscript_DGE_reviewer_response_00_cache_subclasses.R first.
#
# Cell-level negative-binomial mixed model for the IEG panel:
#   counts ~ activity_condition + sex + offset(log(nCount_RNA)) + (1|sample)
# Fixed-effect contrasts vs. SE give log2FC estimates whose SEs account for
# nuclei being nested within animals (each `sample` is one animal).
#
# Only conditions with >= 30 nuclei in this subclass are included (same
# cell_cutoff as GetSubclassContrasts(), used by the DESeq2/MAST/voom
# scripts), so an underpowered condition doesn't get a GLMM fit it can't
# support.

library(Seurat)
library(tidyverse)
library(glmmTMB)
library(broom.mixed)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/reviewer_response"
cache_dir <- file.path(out_dir, "subclass_cache")
dir.create(file.path(out_dir, "glmm"), showWarnings = FALSE, recursive = TRUE)

subclass_list <- LoadSubclassesToUse()
ieg_panel <- LoadGeneList("IEG")

fitGeneGLMM <- function(gene, cell_df, counts_vec) {
  cell_df$count <- counts_vec
  fit <- tryCatch(
    suppressWarnings(
      glmmTMB(count ~ activity_condition + sex + (1 | sample),
              offset = log(nCount_RNA),
              family = nbinom2,
              data = cell_df)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  broom.mixed::tidy(fit, effects = "fixed", conf.int = TRUE) |>
    filter(str_detect(term, "^activity_condition")) |>
    transmute(
      gene = gene,
      contrast = paste0(str_remove(term, "activity_condition"), "_vs_SE"),
      log2FoldChange = estimate / log(2),
      log2FC_CI_low  = conf.low / log(2),
      log2FC_CI_high = conf.high / log(2),
      pvalue = p.value
    )
}

runGeneGLMMs <- function(nuclei_subclass, subclass) {
  df_contrasts <- GetSubclassContrasts(nuclei_subclass, subclass, cell_cutoff = 30) |>
    filter(group2 == "SE")
  valid_conditions <- c("SE", as.character(df_contrasts$group1))
  if (!"SE" %in% valid_conditions || length(valid_conditions) < 2) return(tibble())

  cell_df <- nuclei_subclass@meta.data |>
    as_tibble(rownames = "cell") |>
    filter(activity_condition %in% valid_conditions) |>
    select(cell, sample, activity_condition, sex, nCount_RNA) |>
    mutate(activity_condition = droplevels(
      factor(activity_condition,
             levels = c("SE", "EE30m", "EE6h", "KA30m", "KA6h"))))

  counts <- GetAssayData(nuclei_subclass, assay = "RNA", layer = "counts")
  genes <- intersect(ieg_panel, rownames(counts))

  map_dfr(genes, function(g) {
    res <- fitGeneGLMM(g, cell_df, counts[g, cell_df$cell])
    if (is.null(res)) return(tibble())
    mutate(res, subclass = subclass) |> relocate(subclass)
  })
}

df_glmm_all <- map_dfr(seq_along(subclass_list), function(i) {
  subclass <- subclass_list[i]
  message(glue("[{i}/{length(subclass_list)}] NB-GLMM: {subclass}"))
  subclass_fname <- ShrinkSubclassName(subclass)

  nuclei_subclass <- readRDS(file.path(cache_dir, glue("{subclass_fname}.Rds")))
  res <- runGeneGLMMs(nuclei_subclass, subclass)
  rm(nuclei_subclass)
  gc()
  res
})

# BH correction within each subclass x contrast family
df_glmm_all <- df_glmm_all |>
  group_by(subclass, contrast) |>
  mutate(padj = p.adjust(pvalue, method = "BH")) |>
  ungroup()

write_csv(df_glmm_all, file.path(out_dir, "glmm", "IEG_NB-GLMM_animal-random-effect.csv"))
print(glue("Gene-expression GLMM complete: {out_dir}/glmm/IEG_NB-GLMM_animal-random-effect.csv"))
