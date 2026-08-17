## ----SuppFig_ActiveCellSensitivity
# Reviewer response: sensitivity analysis for the "active cell" definition used in Figure 4.
#
# The main analysis (Fig 4) defines an active cell as one with >= 3 IEGs expressed
# above their 90th-percentile threshold in the SE condition, computed via FindActiveCells().
#
# PART I:  Parameter sweep (consistent with FindActiveCells approach)
#   - Percentile cutoff:  80, 85, 90, 95%
#   - Min. IEGs required: 2, 3, 4
#   Output: gt tables + heatmap of fraction active per supertype × condition × params
#
# PART II: Continuous IEG activity score (mean log-normalised SCT expression across
#          the IEG panel — analogous to AUCell / ssGSEA module scoring)
#   Concordance with the binary call (90th pct, ≥3 IEGs):
#   (a) ROC curve + AUC per condition
#   (b) PR  curve + AUPRC per condition
#   (c) Per-supertype ROC-AUC across conditions
#   (d) Continuous score distribution (active vs. inactive cells) per supertype

library(gt)
library(gtExtras)

source("03-scripts/R/seq_functions.R")

if (!exists("nuclei_ca1")) {
  nuclei <- LoadDataset("Dec2024")
  nuclei_ca1 <- subset(nuclei, subclass_name == "016 CA1-ProS Glut")
  rm(nuclei)  # free up memory; not needed for this analysis
}

supertype_colors <- LoadAllenColors("supertype")
activity_colors  <- LoadActivityColors()

SAVE_PLOTS <- TRUE
save_dir   <- "05-results/SuppFigure_ActiveCellSensitivity/raw_R_plots"

# ── Shared parameters ────────────────────────────────────────────────────────

percentile_cutoffs <- c(0.80, 0.85, 0.90, 0.95)
gene_thresholds    <- c(2L, 3L, 4L)
IEG_symbols        <- LoadGeneList("IEG")
conditions         <- c("SE", "EE30m", "KA30m")

supertype_map <- c(
  "0069 CA1-ProS Glut_1" = "CA1 69-1",
  "0070 CA1-ProS Glut_2" = "CA1 70-2",
  "0071 CA1-ProS Glut_3" = "CA1 71-3",
  "0072 CA1-ProS Glut_4" = "CA1 72-4",
  "0073 CA1-ProS Glut_5" = "CA1 73-5",
  "0074 CA1-ProS Glut_6" = "CA1 74-6"
)
supertype_order     <- unname(supertype_map)
supertype_order_all <- c(supertype_order, "All CA1")

# Supertype colors keyed by short label (consistent with Figure 4D-E)
ca1_supertype_colors <- supertype_colors[names(supertype_map)]
names(ca1_supertype_colors) <- supertype_map
# Color for the combined-subclass row in the heatmap y-axis
ca1_supertype_colors_all <- c(ca1_supertype_colors, "All CA1" = "black")

# ── FindActiveCells: reference call (consistent with Figure 4D-E) ────────────
# gene_threshold = 3, percentile = 90th (hardcoded in FindActiveCells)
if (!exists("FindActiveCells_CA1_outputs")) {
  FindActiveCells_CA1_outputs <- FindActiveCells(
    nuclei_ca1,
    subclass       = "016 CA1-ProS Glut",
    gene_list      = IEG_symbols,
    gene_threshold = 3L
  )
}

df_active      <- FindActiveCells_CA1_outputs$df_active_cells
ref_thresholds <- FindActiveCells_CA1_outputs$activation_thresholds
gene_cols      <- names(ref_thresholds)

# Add supertype info and filter to conditions of interest
supertype_lookup <- nuclei_ca1@meta.data |>
  as_tibble(rownames = "cell") |>
  dplyr::select(cell, supertype_name)

df_active <- df_active |>
  left_join(supertype_lookup, by = "cell") |>
  filter(activity_condition %in% conditions) |>
  mutate(
    activity_condition = factor(activity_condition, levels = conditions),
    supertype_label    = recode(supertype_name, !!!supertype_map)
  )

# ── Helper: recompute per-gene thresholds at a different percentile ───────────
# Thresholds are derived from SE-condition cells in df_active, matching the
# same expression data used by FindActiveCells for the 90th-pct reference.

se_expr_long <- df_active |>
  filter(activity_condition == "SE") |>
  dplyr::select(cell, all_of(gene_cols)) |>
  pivot_longer(cols = -cell, names_to = "gene", values_to = "expression")

compute_thresholds_p <- function(se_long_df, gene_cols_order, p) {
  thresholds <- se_long_df |>
    group_by(gene) |>
    summarise(threshold = quantile(expression, probs = p), .groups = "drop") |>
    arrange(match(gene, gene_cols_order))
  setNames(thresholds$threshold, thresholds$gene)
}

# ═══════════════════════════════════════════════════════════════════════════════
# PART I – Parameter sweep
# ═══════════════════════════════════════════════════════════════════════════════

# Pre-compute n_upregd for each percentile.
# For 90th pct: use num_upregd_genes from FindActiveCells directly (exact match).
# For other percentiles: recompute thresholds from SE cells and count.

n_upregd_per_pct <- list()

for (p in percentile_cutoffs) {
  if (p == 0.90) {
    n_upregd_per_pct[["0.9"]] <- df_active$num_upregd_genes
  } else {
    thresh <- compute_thresholds_p(se_expr_long, gene_cols, p)
    n_upregd_per_pct[[as.character(p)]] <- df_active |>
      rowwise() |>
      mutate(n = sum(c_across(all_of(gene_cols)) > thresh)) |>
      pull(n)
  }
}

# Sweep all 12 parameter combinations
results_list <- list()

for (p in percentile_cutoffs) {
  n_upregd <- n_upregd_per_pct[[as.character(p)]]

  for (gt_n in gene_thresholds) {
    key <- paste0("p", round(p * 100), "_g", gt_n)
    results_list[[key]] <- df_active |>
      ungroup() |>
      mutate(
        active         = n_upregd >= gt_n,
        percentile     = p,
        gene_threshold = gt_n
      ) |>
      dplyr::select(cell, activity_condition, supertype_name, supertype_label,
                    active, percentile, gene_threshold)
  }
}

df_sensitivity <- bind_rows(results_list)

# Summarise ──────────────────────────────────────────────────────────────────

df_summary <- df_sensitivity |>
  group_by(supertype_label, activity_condition, percentile, gene_threshold) |>
  summarise(
    n_cells     = n(),
    n_active    = sum(active),
    frac_active = n_active / n_cells,
    .groups = "drop"
  ) |>
  mutate(
    pct_label = paste0(n_active, " (", round(frac_active * 100), "%)"),
    param_col = paste0("p", round(percentile * 100), "_g", gene_threshold)
  )

# Add a combined row (all CA1 supertypes pooled — no supertype split)
df_summary_all <- df_sensitivity |>
  group_by(activity_condition, percentile, gene_threshold) |>
  summarise(
    n_cells     = n(),
    n_active    = sum(active),
    frac_active = n_active / n_cells,
    .groups = "drop"
  ) |>
  mutate(
    supertype_label = "All CA1",
    pct_label       = paste0(n_active, " (", round(frac_active * 100), "%)"),
    param_col       = paste0("p", round(percentile * 100), "_g", gene_threshold)
  )

df_summary <- bind_rows(df_summary, df_summary_all)

# gt tables ──────────────────────────────────────────────────────────────────

make_sensitivity_table <- function(df, cond, stype_colors) {

  col_order <- paste0(
    rep(paste0("p", c(80, 85, 90, 95)), each = 3),
    "_g", rep(c(2, 3, 4), times = 4)
  )

  df_wide <- df |>
    filter(activity_condition == cond) |>
    dplyr::select(supertype_label, param_col, pct_label) |>
    pivot_wider(names_from = param_col, values_from = pct_label) |>
    left_join(
      df |>
        filter(activity_condition == cond, percentile == 0.90, gene_threshold == 3) |>
        dplyr::select(supertype_label, n_cells),
      by = "supertype_label"
    ) |>
    dplyr::select(supertype_label, n_cells, all_of(col_order)) |>
    mutate(supertype_label = factor(supertype_label, levels = supertype_order_all)) |>
    arrange(supertype_label)

  rev_map     <- setNames(names(supertype_map), supertype_map)
  plot_colors <- c(
    stype_colors[rev_map[supertype_order]],
    "All CA1" = "black"
  )
  names(plot_colors)[seq_along(supertype_order)] <- supertype_order

  tbl <- df_wide |>
    gt(rowname_col = "supertype_label") |>
    tab_header(
      title    = paste0("Active cell counts by parameter set — ", cond),
      subtitle = "N active (% of supertype total) | threshold set from SE-condition IEG expression"
    ) |>
    tab_spanner(label = "80th percentile",           columns = starts_with("p80")) |>
    tab_spanner(label = "85th percentile",           columns = starts_with("p85")) |>
    tab_spanner(label = "90th percentile (Fig. 4)", columns = starts_with("p90")) |>
    tab_spanner(label = "95th percentile",           columns = starts_with("p95")) |>
    cols_label(
      n_cells = "Total N",
      p80_g2 = "≥2 IEGs", p80_g3 = "≥3 IEGs", p80_g4 = "≥4 IEGs",
      p85_g2 = "≥2 IEGs", p85_g3 = "≥3 IEGs", p85_g4 = "≥4 IEGs",
      p90_g2 = "≥2 IEGs", p90_g3 = "≥3 IEGs", p90_g4 = "≥4 IEGs",
      p95_g2 = "≥2 IEGs", p95_g3 = "≥3 IEGs", p95_g4 = "≥4 IEGs"
    ) |>
    tab_style(
      style     = cell_text(weight = "bold"),
      locations = cells_column_spanners()
    ) |>
    tab_style(
      style     = cell_fill(color = "#F0F4FF"),
      locations = cells_body(columns = p90_g3)
    ) |>
    tab_style(
      style     = cell_fill(color = "#F0F4FF"),
      locations = cells_column_labels(columns = p90_g3)
    ) |>
    tab_footnote(
      footnote  = "Highlighted column (90th pct, ≥3 IEGs) is the parameter set used in Figure 4.",
      locations = cells_column_spanners(spanners = "90th percentile (Fig. 4)")
    ) |>
    tab_style(
      style     = list(cell_text(weight = "bold"), cell_fill(color = "#F8F8F8")),
      locations = cells_stub(rows = supertype_label == "All CA1")
    ) |>
    fmt_number(columns = n_cells, use_seps = TRUE, decimals = 0) |>
    opt_table_font(font = "Helvetica") |>
    cols_align(align = "center", columns = everything()) |>
    cols_align(align = "left",   columns = n_cells)

  # data_color() only targets body columns, not the stub. Apply stub text
  # colors row-by-row via tab_style(cells_stub()) instead.
  for (stype in names(plot_colors)) {
    tbl <- tbl |> tab_style(
      style     = cell_text(color = plot_colors[[stype]]),
      locations = cells_stub(rows = supertype_label == stype)
    )
  }
  tbl
}

tables <- setNames(
  lapply(conditions, \(cond) make_sensitivity_table(df_summary, cond, supertype_colors)),
  conditions
)

print(tables$SE)
print(tables$EE30m)
print(tables$KA30m)

# Heatmap: fraction active across all 12 parameter combinations ──────────────

param_levels <- c(
  "80th pct / ≥2 IEGs", "80th pct / ≥3 IEGs", "80th pct / ≥4 IEGs",
  "85th pct / ≥2 IEGs", "85th pct / ≥3 IEGs", "85th pct / ≥4 IEGs",
  "90th pct / ≥2 IEGs", "90th pct / ≥3 IEGs", "90th pct / ≥4 IEGs",
  "95th pct / ≥2 IEGs", "95th pct / ≥3 IEGs", "95th pct / ≥4 IEGs"
)

df_heatmap <- df_summary |>
  filter(activity_condition %in% conditions) |>
  mutate(
    supertype_label    = factor(supertype_label, levels = rev(supertype_order_all)),
    activity_condition = factor(activity_condition, levels = conditions),
    param_label  = paste0(round(percentile * 100), "th pct / ≥", gene_threshold, " IEGs"),
    param_label  = factor(param_label, levels = param_levels),
    is_reference = (percentile == 0.90 & gene_threshold == 3)
  )

p_sensitivity_heatmap <- ggplot(df_heatmap) +
  aes(x = param_label, y = supertype_label, fill = frac_active) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_tile(
    data  = filter(df_heatmap, is_reference),
    color = "black", linewidth = 1.2, fill = NA
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name   = "Fraction\nactive",
    limits = c(0, 1),
    labels = scales::percent
  ) +
  facet_wrap(~activity_condition, nrow = 1) +
  labs(
    x       = NULL, y = NULL,
    title   = "Fraction of active CA1 cells across parameter sets",
    caption = "Bold outline = parameter set used in Figure 4 (90th pct, ≥3 IEGs)"
  ) +
  theme(
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y     = element_text(color = ca1_supertype_colors_all[rev(supertype_order_all)]),
    strip.text      = element_text(face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

print(p_sensitivity_heatmap)

# ═══════════════════════════════════════════════════════════════════════════════
# PART II – Continuous IEG activity score & concordance
# ═══════════════════════════════════════════════════════════════════════════════
# Continuous score: mean log-normalised SCT expression across all IEG panel genes
# per cell. Derived from the same expression matrix used by FindActiveCells, so
# the data are identical. Analogous to AUCell / ssGSEA module scoring.
#
# The reference binary call (active_binary from FindActiveCells: 90th pct, ≥3 IEGs)
# is used as the "ground truth"; the continuous score is the predictor.
# Concordance is quantified via:
#   (a) ROC curve + AUC per condition
#   (b) PR  curve + AUPRC per condition
#   (c) Per-supertype ROC-AUC across conditions
#   (d) Score distributions for active vs. inactive cells per supertype

library(pROC)    # ROC / AUC
library(PRROC)   # precision-recall curves

df_continuous <- df_active |>
  mutate(
    active_ref = active_binary,                         # reference: 90th pct, ≥3 IEGs
    ieg_score  = rowMeans(dplyr::pick(all_of(gene_cols)))  # continuous IEG module score
  )

# (a) Per-condition ROC curves ────────────────────────────────────────────────

roc_df_list <- list()
for (cond in conditions) {
  d       <- filter(df_continuous, activity_condition == cond)
  r       <- pROC::roc(d$active_ref, d$ieg_score, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(r))
  roc_df_list[[cond]] <- tibble(
    fpr             = 1 - r$specificities,
    tpr             = r$sensitivities,
    auc             = auc_val,
    condition       = cond,
    condition_label = paste0(cond, "  (AUC = ", round(auc_val, 3), ")")
  )
}
roc_df <- bind_rows(roc_df_list) |>
  mutate(condition = factor(condition, levels = conditions))

cond_roc_labels <- roc_df |>
  distinct(condition, condition_label) |>
  arrange(match(condition, conditions)) |>
  pull(condition_label)

p_roc <- ggplot(roc_df) +
  aes(x = fpr, y = tpr, color = condition, group = condition) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    values = activity_colors[conditions],
    labels = cond_roc_labels,
    name   = NULL
  ) +
  scale_x_continuous(labels = scales::percent, expand = c(0.01, 0)) +
  scale_y_continuous(labels = scales::percent, expand = c(0.01, 0)) +
  labs(
    x        = "False positive rate",
    y        = "True positive rate",
    title    = "ROC: continuous IEG score vs. binary active call",
    subtitle = "Binary reference: 90th pct threshold, ≥3 IEGs"
  ) +
  theme(legend.position = c(0.65, 0.2))

print(p_roc)

# (b) Per-condition PR curves ─────────────────────────────────────────────────

pr_df_list <- list()
for (cond in conditions) {
  d  <- filter(df_continuous, activity_condition == cond)
  pr <- PRROC::pr.curve(
    scores.class0 = d$ieg_score[d$active_ref],
    scores.class1 = d$ieg_score[!d$active_ref],
    curve = TRUE
  )
  pr_df_list[[cond]] <- as_tibble(pr$curve) |>
    setNames(c("recall", "precision", "threshold")) |>
    mutate(
      condition       = cond,
      auprc           = pr$auc.integral,
      condition_label = paste0(cond, "  (AUPRC = ", round(pr$auc.integral, 3), ")")
    )
}
pr_df <- bind_rows(pr_df_list) |>
  mutate(condition = factor(condition, levels = conditions))

cond_pr_labels <- pr_df |>
  distinct(condition, condition_label) |>
  arrange(match(condition, conditions)) |>
  pull(condition_label)

p_pr <- ggplot(pr_df) +
  aes(x = recall, y = precision, color = condition, group = condition) +
  geom_line(linewidth = 0.9) +
  scale_color_manual(
    values = activity_colors[conditions],
    labels = cond_pr_labels,
    name   = NULL
  ) +
  scale_x_continuous(labels = scales::percent, expand = c(0.01, 0)) +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1), expand = c(0.01, 0)) +
  labs(
    x        = "Recall",
    y        = "Precision",
    title    = "PR curve: continuous IEG score vs. binary active call",
    subtitle = "Binary reference: 90th pct threshold, ≥3 IEGs"
  ) +
  theme(legend.position = c(0.55, 0.55))

print(p_pr)

# (c) Per-supertype ROC-AUC across conditions ─────────────────────────────────

roc_stype_list <- list()
for (cond in conditions) {
  # Per-supertype
  for (stype in supertype_order) {
    d <- filter(df_continuous, activity_condition == cond, supertype_label == stype)
    if (nrow(d) < 10 || sum(d$active_ref) < 5 || sum(!d$active_ref) < 5) next
    auc_val <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(d$active_ref, d$ieg_score, quiet = TRUE))),
      error = \(e) NA_real_
    )
    roc_stype_list[[paste0(cond, "_", stype)]] <- tibble(
      condition       = cond,
      supertype_label = stype,
      auc_roc         = auc_val,
      n_cells         = nrow(d)
    )
  }
  # Combined (all supertypes pooled)
  d_all <- filter(df_continuous, activity_condition == cond)
  if (sum(d_all$active_ref) >= 5 && sum(!d_all$active_ref) >= 5) {
    auc_val_all <- tryCatch(
      as.numeric(pROC::auc(pROC::roc(d_all$active_ref, d_all$ieg_score, quiet = TRUE))),
      error = \(e) NA_real_
    )
    roc_stype_list[[paste0(cond, "_All")]] <- tibble(
      condition       = cond,
      supertype_label = "All CA1",
      auc_roc         = auc_val_all,
      n_cells         = nrow(d_all)
    )
  }
}
roc_stype_df <- bind_rows(roc_stype_list) |>
  mutate(
    supertype_label = factor(supertype_label, levels = supertype_order_all),
    condition       = factor(condition, levels = conditions)
  )

p_auc_supertype <- ggplot(roc_stype_df) +
  aes(x = condition, y = auc_roc, color = supertype_label, group = supertype_label,
      linewidth = supertype_label == "All CA1") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "gray60") +
  geom_line(alpha = 0.8) +
  geom_point(size = 2.5) +
  scale_linewidth_manual(values = c("TRUE" = 1.4, "FALSE" = 0.6), guide = "none") +
  scale_color_manual(
    values = ca1_supertype_colors_all[supertype_order_all],
    name   = NULL
  ) +
  scale_y_continuous(
    limits = c(0.4, 1.0),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    x     = NULL,
    y     = "ROC-AUC",
    title = "Per-supertype ROC-AUC: continuous IEG score vs. binary call"
  ) +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    legend.position = "right"
  )

print(p_auc_supertype)

# (d) Continuous score distribution by binary call ────────────────────────────

df_continuous_with_all <- bind_rows(
  df_continuous,
  df_continuous |> mutate(supertype_label = "All CA1")
) |>
  mutate(supertype_label = factor(supertype_label, levels = supertype_order_all))

p_score_dist <- ggplot(df_continuous_with_all) +
  aes(x = activity_condition, y = ieg_score,
      fill = active_ref, color = active_ref) +
  geom_violin(alpha = 0.5, scale = "width",
              position = position_dodge(0.8)) +
  geom_boxplot(width = 0.2, alpha = 0.9, outlier.shape = NA,
              #  position = position_dodge(0.8),
               color = "black") +
  scale_fill_manual(
    values = c("FALSE" = "gray70", "TRUE" = "#D81B60"),
    labels = c("FALSE" = "Inactive", "TRUE" = "Active"),
    name   = "Binary call\n(90th pct, ≥3 IEGs)"
  ) +
  scale_color_manual(
    values = c("FALSE" = "gray70", "TRUE" = "#D81B60"),
    guide  = "none"
  ) +
  labs(
    x     = NULL,
    y     = "Mean IEG expression (SCT log-normalised)",
    title = "Continuous IEG activity score by binary call group"
  ) +
  facet_wrap(~supertype_label, nrow = 1, scales = "free_x") +
  theme(
    axis.text.x     = element_text(angle = 30, hjust = 1),
    strip.text      = element_text(face = "bold"),
    legend.position = "right"
  )

print(p_score_dist)

# ═══════════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════════

if (SAVE_PLOTS) {
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  for (cond in conditions) {
    gtsave(tables[[cond]],
           filename = file.path(save_dir, paste0("sensitivity_table_", cond, ".png")))
    gtsave(tables[[cond]],
           filename = file.path(save_dir, paste0("sensitivity_table_", cond, ".pdf")))
  }

  ggsave(plot = p_sensitivity_heatmap,
         path = save_dir, filename = "sensitivity_heatmap.png",
         width = 14, height = 4, dpi = 300)
  ggsave(plot = p_sensitivity_heatmap,
         path = save_dir, filename = "sensitivity_heatmap.svg",
         width = 14, height = 4)

  ggsave(plot = p_roc,
         path = save_dir, filename = "IEG_score_ROC.png",
         width = 5, height = 5, dpi = 300)
  ggsave(plot = p_roc,
         path = save_dir, filename = "IEG_score_ROC.svg",
         width = 5, height = 5)

  ggsave(plot = p_pr,
         path = save_dir, filename = "IEG_score_PR.png",
         width = 5, height = 5, dpi = 300)
  ggsave(plot = p_pr,
         path = save_dir, filename = "IEG_score_PR.svg",
         width = 5, height = 5)

  ggsave(plot = p_auc_supertype,
         path = save_dir, filename = "IEG_score_AUC_supertype.png",
         width = 6, height = 4, dpi = 300)
  ggsave(plot = p_auc_supertype,
         path = save_dir, filename = "IEG_score_AUC_supertype.svg",
         width = 6, height = 4)

  ggsave(plot = p_score_dist,
         path = save_dir, filename = "IEG_score_dist.png",
         width = 12, height = 6, dpi = 300)
  ggsave(plot = p_score_dist,
         path = save_dir, filename = "IEG_score_dist.svg",
         width = 12, height = 6)

  print("Sensitivity analyses saved.")
} else {
  print("Rendered without saving.")
}

## ----