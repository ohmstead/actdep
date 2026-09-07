## ----Figure4_SuppFig2
# Reviewer response: concordance between the binary "active cell" definition used in
# Figure 4 and a continuous IEG activity score.
#
# Continuous score: WeightedActivityScore() — each IEG's log-normalised SCT
# expression is z-scored against the mean/sd of that gene in the scaling condition
# (EE30m), then averaged with equal weights across the panel. This corrects for the
# large differences in baseline transcript abundance between IEGs (e.g. Fos vs
# Nr4a1) that a plain row mean would let dominate the score. Derived from the same
# expression matrix used by FindActiveCells, so the underlying data are identical.
#
# The reference binary call (active_binary from FindActiveCells: 90th pct, ≥3 IEGs)
# is used as the "ground truth"; the continuous score is the predictor.
# Concordance is quantified via:
#   (a) ROC curve + AUC per condition
#   (b) PR  curve + AUPRC per condition
#   (c) Per-supertype ROC-AUC across conditions
#   (d) Score distributions for active vs. inactive cells per supertype
#
# The companion parameter sweep (percentile cutoff × min. IEGs required) lives in
# Figure4_SuppFig1.R.

library(pROC)    # ROC / AUC
library(PRROC)   # precision-recall curves

source("03-scripts/R/seq_functions.R")

if (!exists("nuclei_ca1")) {
  # nuclei <- LoadDataset("Dec2024")
  # nuclei_ca1 <- subset(nuclei, subclass_name == "016 CA1-ProS Glut")
  # rm(nuclei)  # free up memory; not needed for this analysis
  nuclei_ca1 <- readRDS("04-analysis/Seurats/nuclei_ca1.rds")
}

supertype_colors <- LoadAllenColors("supertype")
activity_colors  <- LoadActivityColors()

SAVE_PLOTS <- FALSE
save_dir   <- "05-results/Figure4_SuppFig2/raw_R_plots"

# ── Shared parameters ────────────────────────────────────────────────────────

IEG_symbols <- LoadGeneList("IEG")
conditions  <- c("SE", "EE30m", "KA30m")

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
# Color for the combined-subclass series
ca1_supertype_colors_all <- c(ca1_supertype_colors, "All CA1" = "black")

# ── FindActiveCells: reference call (consistent with Figure 4D-E) ────────────
# gene_threshold = 3, percentile = 90th (hardcoded in FindActiveCells)
FindActiveCells_CA1_outputs <- FindActiveCells(
  nuclei_ca1,
  subclass       = "016 CA1-ProS Glut",
  gene_list      = IEG_symbols,
  gene_threshold = 3L
)

df_active <- FindActiveCells_CA1_outputs$df_active_cells

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

# ═══════════════════════════════════════════════════════════════════════════════
# Continuous IEG activity score & concordance
# ═══════════════════════════════════════════════════════════════════════════════

# Continuous, baseline-corrected IEG score (z-scored per gene against SCALE_CONDITION)
WeightedActivityScore_CA1_outputs <- WeightedActivityScore(
  nuclei_ca1,
  subclass        = "016 CA1-ProS Glut",
  gene_list       = IEG_symbols,
  scale_condition = "SE"
)

df_weighted_score <- WeightedActivityScore_CA1_outputs$df_weighted_score

df_continuous <- df_active |>
  mutate(active_ref = active_binary) |>                 # reference: 90th pct, ≥3 IEGs
  left_join(
    df_weighted_score |> dplyr::select(cell, weighted_ieg_score),
    by = "cell"
  ) |>
  rename(ieg_score = weighted_ieg_score)

# (a) Per-condition ROC curves ────────────────────────────────────────────────

roc_df_list <- list()
for (cond in conditions) {
  d       <- filter(df_continuous, activity_condition == cond)
  r       <- pROC::roc(d$active_ref, d$ieg_score, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(r))
  # pROC returns the curve in descending fpr/tpr order. Sort ascending so the
  # points are in curve order for geom_path(); without this, a condition with
  # (near-)perfect separation collapses into one long fpr == 0 tie block and the
  # path jumps from (0, 0) to the first fpr > 0 point as a spurious diagonal.
  roc_df_list[[cond]] <- tibble(
    fpr             = 1 - r$specificities,
    tpr             = r$sensitivities,
    auc             = auc_val,
    condition       = cond,
    condition_label = paste0(cond, "  (AUC = ", sprintf("%.3f", auc_val), ")")
  ) |>
    arrange(fpr, tpr)
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
  geom_path(linewidth = 0.9) +
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
    title    = "ROC: continuous IEG score vs. binary assignment",
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
      condition_label = paste0(cond, "  (AUPRC = ", sprintf("%.3f", pr$auc.integral), ")")
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
  geom_path(linewidth = 0.9) +
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
    title    = "PR curve: continuous IEG score vs. binary assignment",
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
    title = "Per-supertype ROC-AUC: weighted IEG score vs. binary call"
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
               position = position_dodge(0.8),
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
    y     = "Weighted IEG activity score (mean z-score)",
    title = "Weighted IEG activity score by binary call group"
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

  ggsave(plot = p_roc,
         path = save_dir, filename = "IEG_score_ROC.png",
         width = 5.5, height = 5, dpi = 300)
  ggsave(plot = p_roc,
         path = save_dir, filename = "IEG_score_ROC.svg",
         width = 5.5, height = 5)

  ggsave(plot = p_pr,
         path = save_dir, filename = "IEG_score_PR.png",
         width = 5.5, height = 5, dpi = 300)
  ggsave(plot = p_pr,
         path = save_dir, filename = "IEG_score_PR.svg",
         width = 5.5, height = 5)

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

  print("Concordance analyses saved.")
} else {
  print("Rendered without saving.")
}
