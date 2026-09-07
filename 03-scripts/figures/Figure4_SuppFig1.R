## ----Figure4_SuppFig1
# Reviewer response: sensitivity analysis for the "active cell" definition used in Figure 4.
#
# The main analysis (Fig 4) defines an active cell as one with >= 3 IEGs expressed
# above their 90th-percentile threshold in the SE condition, computed via FindActiveCells().
#
# This script: parameter sweep (consistent with the FindActiveCells approach)
#   - Percentile cutoff:  80, 85, 90, 95%
#   - Min. IEGs required: 2, 3, 4
#   Output: gt tables + heatmap of fraction active per supertype × condition × params
#
# The companion concordance analysis (continuous IEG activity score vs. the binary
# call: ROC/PR curves and score distributions) lives in Figure4_SuppFig2.R.

library(gt)
library(gtExtras)

source("03-scripts/R/seq_functions.R")

if (!exists("nuclei_ca1")) {
  # nuclei <- LoadDataset("Dec2024")
  # nuclei_ca1 <- subset(nuclei, subclass_name == "016 CA1-ProS Glut")
  # rm(nuclei)  # free up memory; not needed for this analysis
  nuclei_ca1 <- readRDS("04-analysis/Seurats/nuclei_ca1.rds")
}

supertype_colors <- LoadAllenColors("supertype")

SAVE_PLOTS <- FALSE

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
# Parameter sweep
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
    # title   = "Fraction of active CA1 cells across parameter sets",
    x = NULL, 
    y = NULL
  ) +
  theme(
    # plot.title      = element_text(size = 12),
    strip.text.x.top = element_text(size = 12, face = "bold"),
    axis.text.x     = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y     = element_text(color = ca1_supertype_colors_all[rev(supertype_order_all)], size = 12),
    strip.text      = element_text(face = "bold"),
    panel.grid      = element_blank(),
    legend.position = "right"
  )

print(p_sensitivity_heatmap)

# ═══════════════════════════════════════════════════════════════════════════════
# Save
# ═══════════════════════════════════════════════════════════════════════════════

if (SAVE_PLOTS) {
  save_dir   <- "05-results/SuppItem2/raw_R_plots"
  dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)

  for (cond in conditions) {
    gtsave(tables[[cond]],
           filename = file.path(save_dir, paste0("sensitivity_table_", cond, ".png")))
    gtsave(tables[[cond]],
           filename = file.path(save_dir, paste0("sensitivity_table_", cond, ".pdf")))
  }
  
  save_dir   <- "05-results/Figure4_SuppFig1/raw_R_plots"
  ggsave(plot = p_sensitivity_heatmap,
         path = save_dir, filename = "sensitivity_heatmap.png",
         width = 14, height = 4, dpi = 300)
  ggsave(plot = p_sensitivity_heatmap,
         path = save_dir, filename = "sensitivity_heatmap.svg",
         width = 14, height = 4)

  print("Sensitivity analyses saved.")
} else {
  print("Rendered without saving.")
}
