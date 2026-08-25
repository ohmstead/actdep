# Reviewer response follow-up: rank-rank hypergeometric overlap (RRHO)
# comparing DESeq2 pseudobulk vs. the NB-GLMM (animal random effect) gene
# rankings for one subclass x contrast, styled after Figure 1A/1D of
# Piron et al. 2024, Life Science Alliance (RedRibbon paper,
# https://www.life-science-alliance.org/content/7/2/e202302203):
#   - a "real data" map showing a diagonal signal from down- to
#     up-regulation when the two rankings agree (like their Fig 1A, built
#     from two identical lists)
#   - a negative-control map with one list's gene order scrambled, which
#     should show no signal (like their Fig 1D, built from random lists)
#
# Subclass/contrast: 016 CA1-ProS Glut, EE30m vs. SE.
#
# The NB-GLMM in manuscript_DGE_reviewer_response_03_glmm_gene_expression.R
# was deliberately restricted to the 15-gene IEG panel -- fitting it
# genome-wide (glmmTMB, per gene, ~20,000 genes passing the 1%-expression
# filter in CA1) benchmarks at ~2.8 sec/gene, i.e. ~15 hours serial for this
# one subclass x contrast alone. Per user direction, this script instead
# fits the top 2000 genes by mean expression (not the full transcriptome)
# and parallelizes the fits across cores to make that tractable (~2000
# genes / 8 cores * 2.8 sec =~ 12 min).
#
# RRHO method: genes are ranked within each list by a signed statistic
# (sign(log2FC) * -log10(p)), from most down-regulated to most
# up-regulated. At a grid of rank-cutoff pairs (i, j), the number of genes
# in common between the top-i of list 1 and top-j of list 2 is tested
# against a hypergeometric null; the resulting -log10(p) grid is the RRHO
# map. This is the same test Plaisier et al. 2010 introduced and Piron et
# al. 2024 describes in their Methods.

library(Seurat)
library(tidyverse)
library(glmmTMB)
library(parallel)

source("03-scripts/R/seq_functions.R")

out_dir <- "04-analysis/reviewer_response"
cache_dir <- file.path(out_dir, "subclass_cache")
deseq_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
rrho_dir <- file.path(out_dir, "rrho")
dir.create(rrho_dir, showWarnings = FALSE, recursive = TRUE)

subclass <- "016 CA1-ProS Glut"
subclass_fname <- ShrinkSubclassName(subclass)
group1 <- "EE30m"
group2 <- "SE"
n_top_genes <- 2000
n_cores <- if (.Platform$OS.type == "windows") max(1, min(4, detectCores() - 1)) else max(1, min(8, detectCores() - 1))


# ============================================================================ #
# genome-wide (top-N by expression) NB-GLMM for this one subclass x contrast ----
# ============================================================================ #
ca1 <- readRDS(file.path(cache_dir, glue("{subclass_fname}.Rds")))
ca1 <- subset(ca1, activity_condition %in% c(group1, group2))

mat <- ca1[["SCT"]]@data
expressed_genes <- rownames(mat)[rowMeans(mat > 0) > 0.01]
mean_expr <- rowMeans(mat[expressed_genes, , drop = FALSE])
top_genes <- names(sort(mean_expr, decreasing = TRUE))[seq_len(n_top_genes)]
rm(mat, mean_expr)

cell_df <- ca1@meta.data |>
  as_tibble(rownames = "cell") |>
  select(cell, sample, activity_condition, nCount_RNA) |>
  mutate(activity_condition = droplevels(factor(activity_condition, levels = c(group2, group1))))

counts <- GetAssayData(ca1, assay = "RNA", layer = "counts")
counts_sub <- as.matrix(counts[top_genes, cell_df$cell, drop = FALSE])
rm(ca1, counts)
gc()

fitGeneGLMM <- function(gene) {
  df <- cell_df
  df$count <- counts_sub[gene, df$cell]
  fit <- tryCatch(
    suppressWarnings(
      glmmTMB(count ~ activity_condition + (1 | sample),
              offset = log(nCount_RNA), family = nbinom2, data = df)),
    error = function(e) NULL
  )
  if (is.null(fit)) return(tibble(gene = gene, log2FoldChange = NA_real_, pvalue = NA_real_))

  co <- summary(fit)$coefficients$cond
  term <- grep("^activity_condition", rownames(co), value = TRUE)
  if (length(term) == 0) return(tibble(gene = gene, log2FoldChange = NA_real_, pvalue = NA_real_))

  tibble(gene = gene,
         log2FoldChange = co[term, "Estimate"] / log(2),
         pvalue = co[term, "Pr(>|z|)"])
}

message(glue("Fitting genome-wide (top {n_top_genes}) NB-GLMM on {n_cores} cores: {subclass} {group1} vs {group2}"))
t0 <- Sys.time()

batch_size <- 100
gene_batches <- split(top_genes, ceiling(seq_along(top_genes) / batch_size))

cl <- makeCluster(n_cores)
clusterEvalQ(cl, { library(glmmTMB); library(tibble) })
clusterExport(cl, c("cell_df", "counts_sub"))

n_done <- 0
glmm_results <- list()
for (batch in gene_batches) {
  glmm_results <- c(glmm_results, parLapply(cl, batch, fitGeneGLMM))
  n_done <- n_done + length(batch)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  eta_min <- round((elapsed / n_done) * (n_top_genes - n_done) / 60, 1)
  message(glue("  [{n_done}/{n_top_genes}] genes fit, {round(elapsed / 60, 1)} min elapsed, ~{eta_min} min remaining"))
}
stopCluster(cl)

message(glue("GLMM fits complete in {round(difftime(Sys.time(), t0, units = 'mins'), 1)} min"))

df_glmm <- bind_rows(glmm_results) |> drop_na(log2FoldChange, pvalue)
write_csv(df_glmm, file.path(rrho_dir, glue("{subclass_fname}__{group1}_vs_{group2}__NB-GLMM_top{n_top_genes}genes.csv")))

rm(counts_sub, cell_df, glmm_results)
gc()


# ============================================================================ #
# signed rankings, DESeq2 vs. NB-GLMM, shared gene set ----
# ============================================================================ #
# computeRRHO()/plotRRHO() are shared helpers, sourced from seq_functions.R above.
deseq <- read_csv(file.path(deseq_dir, glue("{subclass_fname}__{group1}_vs_{group2}.csv")),
                  show_col_types = FALSE) |>
  filter(gene %in% top_genes) |>
  drop_na(log2FoldChange.shrink, pvalue) |>
  transmute(gene, signed_stat = sign(log2FoldChange.shrink) * -log10(pmax(pvalue, 1e-300)))

glmm <- df_glmm |>
  transmute(gene, signed_stat = sign(log2FoldChange) * -log10(pmax(pvalue, 1e-300)))

# --- real-data RRHO map ----
df_rrho_real <- computeRRHO(deseq, glmm)
write_csv(df_rrho_real, file.path(rrho_dir, glue("{subclass_fname}__{group1}_vs_{group2}__RRHO_real.csv")))

p_real <- plotRRHO(df_rrho_real,
                   glue("RRHO: DESeq2 vs. NB-GLMM ({subclass}, {group1} vs {group2})"))
ggsave(file.path(rrho_dir, "RRHO_real.pdf"), p_real, width = 6, height = 5.5)

# --- negative control: scramble the GLMM gene order ----
set.seed(17)
glmm_scrambled <- glmm |>
  mutate(signed_stat = sample(signed_stat))

df_rrho_scrambled <- computeRRHO(deseq, glmm_scrambled)
write_csv(df_rrho_scrambled, file.path(rrho_dir, glue("{subclass_fname}__{group1}_vs_{group2}__RRHO_scrambled.csv")))

p_scrambled <- plotRRHO(df_rrho_scrambled,
                        glue("RRHO negative control: DESeq2 vs. scrambled NB-GLMM ({subclass}, {group1} vs {group2})"))
ggsave(file.path(rrho_dir, "RRHO_scrambled_negative_control.pdf"), p_scrambled, width = 6, height = 5.5)

print(glue("RRHO analysis complete. Outputs in {rrho_dir}/"))
