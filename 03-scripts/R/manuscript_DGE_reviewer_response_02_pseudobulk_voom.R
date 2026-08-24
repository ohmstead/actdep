# Reviewer response (b): complementary pseudobulk DGE framework (limma-voom),
# plus cross-method concordance against the existing DESeq2 pseudobulk results.
# Run manuscript_DGE_reviewer_response_00_cache_subclasses.R first.
#
# Same unit of analysis as the DESeq2 pseudobulk models (animal-level counts),
# with sex included as a covariate. All contrasts from GetSubclassContrasts()
# are tested. voom results include 95% CIs on logFC for the forest plots.
#
# NOTE: an earlier version of this script also ran edgeR-QLF as a second
# complementary framework. Dropped in favor of limma-voom alone: voom's
# topTable(confint=TRUE) reports a logFC standard error/CI directly, which
# is what the reviewer's forest-plot request needs; edgeR-QLF's glmQLFTest
# doesn't expose that without extra derivation. voom's mean-variance
# weighting also suits the pseudobulk library-size variation across
# subclasses (n=2-15 animals per condition) well.

library(Seurat)
library(tidyverse)
library(edgeR)
library(limma)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/reviewer_response"
cache_dir <- file.path(out_dir, "subclass_cache")
deseq_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
dir.create(file.path(out_dir, "voom"), showWarnings = FALSE, recursive = TRUE)

subclass_list <- LoadSubclassesToUse()
log2FC_threshold <- 0.585   # 1.5 fold-change


buildPseudobulk <- function(nuclei_subclass) {
  # only test genes expressed in >= 1% of nuclei (matches DESeq2 script)
  mat <- nuclei_subclass[["SCT"]]@data
  gene_list <- rownames(mat)[rowMeans(mat > 0) > 0.01]

  counts <- AggregateExpression(nuclei_subclass,
                                group.by = "sample",
                                features = gene_list,
                                return.seurat = FALSE)$RNA

  sample_meta <- nuclei_subclass@meta.data |>
    as_tibble() |>
    distinct(sample, activity_condition, sex) |>
    mutate(sample = str_replace_all(sample, "_", "-")) |>
    column_to_rownames("sample")
  sample_meta <- sample_meta[colnames(counts), ]

  list(counts = as.matrix(counts), meta = sample_meta)
}

runVoom <- function(pb, subclass_fname, df_contrasts) {
  condition <- droplevels(factor(pb$meta$activity_condition))
  sex <- factor(pb$meta$sex)
  design <- model.matrix(~ 0 + condition + sex)
  colnames(design) <- sub("^condition", "", colnames(design))

  y <- DGEList(counts = pb$counts, group = condition)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  v <- voom(y, design)
  fit_voom <- lmFit(v, design)

  for (k in seq_len(nrow(df_contrasts))) {
    group1 <- as.character(df_contrasts$group1[k])
    group2 <- as.character(df_contrasts$group2[k])
    if (!all(c(group1, group2) %in% colnames(design))) next
    contrast_vec <- makeContrasts(
      contrasts = glue("{group1} - {group2}"), levels = design)

    # limma-voom (with 95% CI on logFC)
    fit_c <- contrasts.fit(fit_voom, contrast_vec) |> eBayes(robust = TRUE)
    res_voom <- topTable(fit_c, number = Inf, confint = 0.95) |>
      rownames_to_column("gene") |>
      as_tibble() |>
      rename(log2FoldChange = logFC, log2FC_CI_low = CI.L, log2FC_CI_high = CI.R,
             pvalue = P.Value, padj = adj.P.Val) |>
      mutate(classification = case_when(
        log2FoldChange >  log2FC_threshold & padj < 0.05 ~ "upregulated",
        log2FoldChange < -log2FC_threshold & padj < 0.05 ~ "downregulated",
        TRUE ~ "no_change"
      ))
    write_csv(res_voom,
              glue("{out_dir}/voom/{subclass_fname}__{group1}_vs_{group2}.csv"))
  }
}

for (i in seq_along(subclass_list)) {
  subclass <- subclass_list[i]
  message(glue("[{i}/{length(subclass_list)}] pseudobulk voom: {subclass}"))
  subclass_fname <- ShrinkSubclassName(subclass)

  nuclei_subclass <- readRDS(file.path(cache_dir, glue("{subclass_fname}.Rds")))
  df_contrasts <- GetSubclassContrasts(nuclei_subclass, subclass, cell_cutoff = 30)
  pseudo <- buildPseudobulk(nuclei_subclass)

  rm(nuclei_subclass)
  gc()

  runVoom(pseudo, subclass_fname, df_contrasts)
  rm(pseudo)
}


# --- cross-method concordance with existing DESeq2 results ----
loadDEGSet <- function(path) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) |>
    filter(classification != "no_change") |>
    pull(gene)
}

df_concordance <- map_dfr(subclass_list, function(subclass) {
  subclass_fname <- ShrinkSubclassName(subclass)
  contrast_files <- list.files(file.path(out_dir, "voom"),
                               pattern = glue("^{subclass_fname}__"),
                               full.names = FALSE)
  map_dfr(contrast_files, function(f) {
    degs <- list(
      DESeq2 = loadDEGSet(file.path(deseq_dir, f)),
      voom   = loadDEGSet(file.path(out_dir, "voom", f))
    )
    degs <- degs[!map_lgl(degs, is.null)]
    jaccard <- function(a, b) {
      if (length(union(a, b)) == 0) return(NA_real_)
      length(intersect(a, b)) / length(union(a, b))
    }
    tibble(
      subclass = subclass,
      contrast = str_remove(f, glue("^{subclass_fname}__")) |> str_remove(".csv$"),
      n_DEG_DESeq2 = length(degs$DESeq2),
      n_DEG_voom   = length(degs$voom),
      n_DEG_both_methods = length(reduce(degs, intersect)),
      jaccard_DESeq2_voom = jaccard(degs$DESeq2, degs$voom)
    )
  })
})

write_csv(df_concordance, file.path(out_dir, "0_method_concordance.csv"))
print(df_concordance, n = Inf)
print(glue("Pseudobulk voom + concordance complete. Outputs in {out_dir}/voom/"))
