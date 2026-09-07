## ----Figure5_SuppFig2
# Reviewer response: ZT differential expression in SE re-run with limma-voom
# instead of DESeq2, and the Figure 5B-E heatmaps recreated from that fit.
#
# Everything about the design is held constant against Figure5B-E_heatmaps.R --
# same 4 subclasses (CA1, DG, astrocytes, oligodendrocytes), same SE-only
# nuclei, same animal-level pseudobulk, same 1%-of-nuclei gene filter, same
# 6 pairwise ZT contrasts, same |log2FC| >= 0.585 & FDR < 0.05 DEG call. The
# only thing that changes is the inference framework:
#
#   DESeq2 (Fig 5)      : ~ ZT, negative-binomial GLM, median-of-ratios counts
#   limma-voom (here)   : ~ 0 + ZT, TMM + voom log2-CPM, precision-weighted LM
#
# As in Figure 3 Supp Fig 1, the *visualized* values come from the same
# normalization as the model that called the DEGs: the heatmap body and the tau
# annotation are both derived from the voom log2-CPM matrix (v$E), never from
# DESeq2 normalized counts.
#   - heatmap body: voom log2-CPM, averaged within ZT, z-scored per gene
#   - tau         : voom log2-CPM back-transformed to linear CPM (2^E), so that
#                   tau() operates on a positive abundance scale, as it does on
#                   DESeq2 normalized counts in Figure 5

source("03-scripts/R/seq_functions.R")

library(Seurat)
library(purrr)
library(limma)
library(edgeR)
library(ComplexHeatmap)
library(circlize)
library(seriation)

if (!exists("SAVE_PLOTS")) SAVE_PLOTS <- FALSE

deg_save_path <- "04-analysis/DEGs/Dec2024_ZT_SE_pseudobulk_voom"
deseq_deg_path <- "04-analysis/DEGs/Dec2024_ZT_SE_pseudobulk"
dir.create(deg_save_path, showWarnings = FALSE, recursive = TRUE)

zt_levels <- c("ZT0", "ZT4", "ZT12", "ZT16")
thresh <- 0.585   # 1.5 fold-change, matches Figure 5

zt_contrasts <- tibble::tribble(
  ~group1,  ~group2,
  "ZT4",    "ZT0",
  "ZT12",   "ZT0",
  "ZT16",   "ZT0",
  "ZT12",   "ZT4",
  "ZT16",   "ZT4",
  "ZT16",   "ZT12"
)

if (!exists('nuclei')) {nuclei <- LoadDataset('Dec2024')}

seurat_subsets <- list(
  CA1   = subset(nuclei, subclass_name == '016 CA1-ProS Glut' & activity_condition %in% c('SE')),
  DG    = subset(nuclei, subclass_name == '037 DG Glut'       & activity_condition %in% c('SE')),
  Astro = subset(nuclei, subclass_name == '319 Astro-TE NN'   & activity_condition %in% c('SE')),
  Oligo = subset(nuclei, subclass_name == '327 Oligo NN'      & activity_condition %in% c('SE'))
)


# define fxns ----------------------------------------
buildPseudobulk <- function(nuclei_subclass) {
  # Animal-level pseudobulk counts + per-sample ZT, identical to the input the
  # DESeq2 version of this analysis receives.
  #
  # only test genes expressed in >1% of nuclei (matches the DESeq2 script)
  mat <- nuclei_subclass[["SCT"]]@data
  gene_list <- rownames(mat)[rowMeans(mat > 0) > 0.01]

  counts <- AggregateExpression(nuclei_subclass,
                                group.by = "sample",
                                features = gene_list,
                                return.seurat = FALSE)$RNA

  col_data <- nuclei_subclass@meta.data |>
    as_tibble() |>
    distinct(sample, ZT) |>
    mutate(sample = str_replace_all(sample, "_", "-"),
           ZT = factor(ZT, levels = zt_levels)) |>
    column_to_rownames("sample")
  col_data <- col_data[colnames(counts), , drop = FALSE]
  col_data$sample <- rownames(col_data)
  stopifnot(!any(is.na(col_data$ZT)))

  list(counts = as.matrix(counts), meta = col_data)
}


runVoomZT <- function(pb, subclass) {
  # Fits ~ 0 + ZT with limma-voom and tests every pairwise ZT contrast.
  # Returns both the per-contrast results and the voom object, because v$E is
  # what the heatmap and tau annotation are built from downstream.
  zt <- droplevels(pb$meta$ZT)
  design <- model.matrix(~ 0 + zt)
  colnames(design) <- sub("^zt", "", colnames(design))

  y <- DGEList(counts = pb$counts, group = zt)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")

  v <- voom(y, design, plot = FALSE)
  fit <- lmFit(v, design)

  res_list <- list()
  for (k in seq_len(nrow(zt_contrasts))) {
    group1 <- zt_contrasts$group1[k]
    group2 <- zt_contrasts$group2[k]
    contrast_name <- glue("{group1}_vs_{group2}")
    if (!all(c(group1, group2) %in% colnames(design))) {
      message(glue("Skipping {subclass} {contrast_name}: group missing from design."))
      next
    }

    contrast_vec <- makeContrasts(contrasts = glue("{group1} - {group2}"),
                                  levels = design)
    fit_c <- contrasts.fit(fit, contrast_vec) |> eBayes(robust = TRUE)

    res <- topTable(fit_c, number = Inf, confint = 0.95) |>
      rownames_to_column("gene") |>
      as_tibble() |>
      dplyr::rename(log2FoldChange = logFC, log2FC_CI_low = CI.L,
                    log2FC_CI_high = CI.R, pvalue = P.Value, padj = adj.P.Val) |>
      mutate(classification = case_when(
        log2FoldChange >  thresh & padj < 0.05 ~ "upregulated",
        log2FoldChange < -thresh & padj < 0.05 ~ "downregulated",
        TRUE ~ "no_change"
      )) |>
      arrange(padj)

    print(glue("{subclass} {contrast_name}: {sum(res$classification != 'no_change')} DEGs"))
    write_csv(res, glue("{deg_save_path}/{subclass}_{contrast_name}.csv"))
    res_list[[contrast_name]] <- res
  }

  list(results = res_list, voom = v, meta = pb$meta)
}


# run loop ----------------------------------------
plots <- list()
voom_fits <- list()
for (subclass in names(seurat_subsets)) {
  nuclei_subclass <- seurat_subsets[[subclass]]
  nuclei_subclass$ZT <- factor(nuclei_subclass$ZT, levels = zt_levels)
  nuclei_subclass$sample <- factor(nuclei_subclass$sample)


  # run voom ----------------------------------------
  pb <- buildPseudobulk(nuclei_subclass)
  fit <- runVoomZT(pb, subclass)
  voom_fits[[subclass]] <- fit

  v <- fit$voom
  col_data <- fit$meta


  # heatmap ----------------------------------------
  # assemble all genes and their contrast of origin
  all_genes <- imap_dfr(fit$results, function(res, contrast_name) {
    res |>
      filter(classification != 'no_change') |>
      mutate(contrast = str_replace(contrast_name, "_vs_", " vs "))
  })

  if (nrow(all_genes) == 0) {
    message(glue("Skipping {subclass}: no DE genes identified for heatmap."))
    next
  }
  write_csv(all_genes, glue("{deg_save_path}/{subclass}_all_DEGs.csv"))

  # remove duplicates
  gene_list <- all_genes |>
    distinct(gene) |>
    pull()

  if (length(gene_list) < 5) {
    message(glue("Skipping {subclass} heatmap: only {length(gene_list)} DEG(s)."))
    next
  }

  # sample -> ZT lookup for the voom expression matrix
  sample_zt <- col_data |>
    as_tibble() |>
    select(sample, ZT)

  # get voom-normalized expression (log2-CPM), averaged within ZT, z-scored
  expr_log2 <- v$E[gene_list, , drop = FALSE]
  gene_mat <- expr_log2 |>
    t() |>
    as.data.frame() |>
    rownames_to_column('sample') |>
    pivot_longer(cols = -sample, names_to = 'gene', values_to = 'expression') |>
    left_join(sample_zt, by = 'sample') |>
    group_by(gene, ZT) |>
    dplyr::summarize(mean_expression = mean(expression), .groups = 'drop') |>
    group_by(gene) |>
    mutate(norm_expression = as.numeric(scale(mean_expression))) |>   # Z-SCORE EXPRESSION
    ungroup() |>
    select(gene, ZT, norm_expression) |>
    pivot_wider(names_from = ZT, values_from = norm_expression) |>
    column_to_rownames('gene') |>
    as.matrix()
  gene_mat <- gene_mat[, intersect(zt_levels, colnames(gene_mat)), drop = FALSE]


  # get tau ----------------------------------------
  # tau() expects a positive abundance scale (Fig 5 feeds it DESeq2 normalized
  # counts), so the voom log2-CPM values are back-transformed to linear CPM.
  df_expression_zt <- (2 ^ expr_log2) |>
    as.data.frame() |>
    rownames_to_column('gene') |>
    pivot_longer(cols = -gene, names_to = 'sample', values_to = 'expression') |>
    left_join(sample_zt, by = 'sample') |>
    group_by(gene, ZT) |>
    dplyr::summarize(expression = mean(expression), .groups = 'drop') |>
    dplyr::mutate(tau = tau(expression), .by = gene) |>
    arrange(dplyr::desc(tau)) |>
    print()

  df_expression_zt |>
    slice_max(expression, by = gene) |>
    slice_max(tau, by = ZT) |>
    arrange(ZT)


  # seriate matrix ----------------------------------------
  o <- seriate(gene_mat, method = 'Heatmap', seriation_method = 'OLO_average') |> get_order(1)
  gene_mat_ordered <- gene_mat[names(o),]
  hc <- hclust(dist(gene_mat), method = "average")
  dend <- as.dendrogram(hc)
  dend <- dendextend::rotate(dend, order = o)


  # make tau anno ----------------------------------------
  tau_thresh <- 0.75
  genes_tyssowski <- c(LoadGeneList('tyssowski'), 'Adcy1', 'Adcy8')
  genes_clock <- LoadGeneList('circadian')

  anno_colors <- setNames(c('red', 'orange', 'black'),
                          c('high_tau', 'clock', 'tyssowski'))
  anno_colors

  # df for tau plot
  df_anno <- df_expression_zt |>
    distinct(gene, tau) |>
    mutate(gene = factor(gene, levels = rownames(gene_mat_ordered))) |>
    arrange(gene) |>
    mutate(color = ifelse(tau > tau_thresh, anno_colors[1],
                          ifelse(gene %in% genes_clock, anno_colors[2],
                                 anno_colors[3]))) |> print()

  # get high tau genes
  high_tau_genes <- df_anno |>
    filter(tau > tau_thresh) |>
    mutate(text_color = anno_colors[1]) |>
    select(gene, text_color) |>
    group_by(gene) |>
    slice_head() |>
    mutate(type = 'high_tau') |>
    print()

  # add any other genes of interest to list
  clock_genes <- tibble(
    gene = genes_clock,
    text_color = anno_colors[2],
    type = 'clock'
  ) |>
    filter(!gene %in% high_tau_genes$gene)
  tyssowski_genes <- tibble(
    gene = genes_tyssowski,
    text_color = anno_colors[3],
    type = 'tyssowski'
  ) |>
    filter(!gene %in% high_tau_genes$gene) |>
    filter(!gene %in% clock_genes$gene)

  df_marked_genes <- rbind(high_tau_genes, clock_genes, tyssowski_genes) |>
    filter(gene %in% rownames(gene_mat_ordered)) |>
    mutate(gene = factor(gene, levels = rownames(gene_mat_ordered))) |>
    arrange(gene) |>
    left_join(read_csv("04-analysis/gene_biotypes.csv", show_col_types = FALSE),
              by = c('gene' = 'gene_name')) |>
    print(n = 31)

  tau_anno <- anno_barplot(
    df_anno$tau,
    bar_width = 0.7,
    width = unit(1, 'in'),
    ylim = c(0,1),
    gp = gpar(fill = df_anno$color, col = df_anno$color),
    which = 'row'
  )

  # anno_mark errors on an empty mark set, which can happen if voom calls a
  # small DEG set with no high-tau/clock/PRG gene in it
  if (nrow(df_marked_genes) > 0) {
    anno_right <- HeatmapAnnotation(
      tau = tau_anno,
      genes = anno_mark(
        at = which(rownames(gene_mat_ordered) %in% df_marked_genes$gene),
        labels = intersect(rownames(gene_mat_ordered), df_marked_genes$gene),
        lines_gp = gpar(col = df_marked_genes$text_color),
        labels_gp = gpar(col = df_marked_genes$text_color, angle = 315),
        side = 'right',
        which = 'row'
      ),
      which = 'row',
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
  } else {
    anno_right <- HeatmapAnnotation(
      tau = tau_anno,
      which = 'row',
      show_annotation_name = FALSE,
      show_legend = FALSE
    )
  }


  # find km clusters ----------------------------------------
  # Seeds/orders carried over from Figure5B-E_heatmaps.R. voom calls a
  # different DEG set than DESeq2 does, so these no longer guarantee the same
  # aesthetic ordering as the main figure and may need re-tuning per subclass.
  subclass_seeds_orders <- tibble(
    subclass_title = c('CA1', 'DG', 'Astro', 'Oligo'),
    seed    = c(17, 1, 9, 5),
    km_order   = list(c(4,1,3,2), c(4,2,3,1), c(1,3,2,4), c(2,4,1,3))
  )
  subclass_specifics <- subclass_seeds_orders |> filter(subclass_title == subclass) |> print()
  seed <- subclass_specifics$seed
  desired_order <- subclass_specifics$km_order[[1]]

  # run km
  n_centers <- min(4, nrow(gene_mat_ordered) - 1)
  set.seed(seed)
  km <- kmeans(gene_mat_ordered, centers = n_centers)
  split <- factor(km$cluster, levels = intersect(desired_order, unique(km$cluster)))


  # plot ----------------------------------------
  width_scalar <- 2
  width_params <- list(
    width = unit(width_scalar, 'in'),
    row_dend_width = unit(width_scalar * 0.075, 'in')
  )
  hm <- Heatmap(
    gene_mat_ordered,
    right_annotation = anno_right,
    row_split = split,
    cluster_rows = T,
    cluster_row_slices = F,  # keeps your factor level order for the slices
    cluster_columns = F,
    show_row_names = F,
    show_column_names = F,
    show_heatmap_legend = F,
    row_title = NULL,
    column_title = glue("{subclass} (limma-voom)"),
    width = width_params$width, # width of the heatmap body
    row_dend_width = width_params$row_dend_width, # width of the row dendrogram
    col = colorRamp2(c(-2,0,2), hcl_palette = "inferno")
  ) |> draw()

  plots <- c(plots, setNames(list(hm), subclass))
}


# DESeq2 vs. voom concordance ----------------------------------------
# Same summary the pseudobulk reviewer-response script reports for the
# activity-condition contrasts, here for the ZT contrasts.
loadDEGSet <- function(path) {
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) |>
    filter(classification != "no_change") |>
    pull(gene)
}

df_concordance <- map_dfr(names(seurat_subsets), function(subclass) {
  map_dfr(seq_len(nrow(zt_contrasts)), function(k) {
    contrast_name <- glue("{zt_contrasts$group1[k]}_vs_{zt_contrasts$group2[k]}")
    fname <- glue("{subclass}_{contrast_name}.csv")
    degs_deseq <- loadDEGSet(file.path(deseq_deg_path, fname))
    degs_voom  <- loadDEGSet(file.path(deg_save_path, fname))
    if (is.null(degs_deseq) || is.null(degs_voom)) return(NULL)

    n_union <- length(union(degs_deseq, degs_voom))
    tibble(
      subclass = subclass,
      contrast = contrast_name,
      n_DEG_DESeq2 = length(degs_deseq),
      n_DEG_voom = length(degs_voom),
      n_DEG_both_methods = length(intersect(degs_deseq, degs_voom)),
      jaccard_DESeq2_voom = if (n_union == 0) NA_real_ else
        length(intersect(degs_deseq, degs_voom)) / n_union
    )
  })
})

write_csv(df_concordance, glue("{deg_save_path}/0_method_concordance_ZT.csv"))
print(df_concordance, n = Inf)


# save ----------------------------------------
if (SAVE_PLOTS) {
  save_path <- "05-results/Figure5_SuppFig2/raw_R_plots"
  dir.create(save_path, showWarnings = FALSE, recursive = TRUE)

  for (subclass in names(plots)) {
    subclass_fname <- ShrinkSubclassName(subclass)
    # png
    png(glue("{save_path}/SE__ZT_heatmap_voom_{subclass_fname}.png"),
        width = 6, height = 10, units = "in", res = 900)
    draw(plots[[subclass]]); dev.off()

    # svg
    svgsave(plot = plots[[subclass]],
            path = save_path,
            filename = glue("SE_ZT_heatmap_voom_{subclass_fname}"),
            width = 6, height = 10)
  }
}
## ----
