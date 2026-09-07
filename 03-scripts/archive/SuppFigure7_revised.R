## ----SuppFig7  —  STEP 1: detection of neurogenesis markers along the stage axis
## ============================================================================
## Manuscript : Olmstead, King & Bloodgood 2025 (eLife) — ACT-DEPP snRNA-seq
##              hippocampus atlas.
## Addresses  : Reviewer #1, recommendation 3 (neurogenesis-related signatures).
##
## DESIGN (agreed in planning):
##   * Genes chosen A PRIORI from the canonical adult-hippocampal-neurogenesis
##     staging scheme of Nicola, Hippenmeyer & Kempermann (2015), Front.
##     Neuroanat. 9:53, Fig 1A. NOT derived from any pseudobulk DEG csv. Each
##     gene is annotated with the stage it marks and a primary reference.
##   * Primary readout = DETECTION along the stage axis (type-1 RGL -> type-2a
##     IPC -> type-2b/3 immature -> mature GC). The point is to show, honestly,
##     which end of the lineage the adult P30 snRNA-seq atlas captures
##     (expected: immature/post-mitotic end yes; rare quiescent/amplifying
##     progenitor end sparse-to-absent), not to assert a trajectory.
##   * Detection quantified from the per-cell RNA assay, SE baseline only:
##       det_frac = fraction of nuclei with raw count > 0
##       mean_exp = mean log-normalized expression across all nuclei in subclass
##     (NOTE: the rest of this codebase computes per-cell expression from the SCT
##      'data' layer — see QuickPercentExpression / GetCorrData. We use the RNA
##      assay here by explicit decision, because the detection claim should rest
##      on raw capture, not SCT-corrected values. ASSAY is a switch below if you
##      want to mirror the codebase default instead.)
##
## SCOPE: STEP 1 only. Steps 2 (eligibility gate), 3 (modulation across
##        condition/ZT), and 4 (narrative + SUPPLEMENT.docx) are added later.
##
## RUN: from project root in the actdep Rproj, as with the other figure scripts.
##      Requires the Dec2024 Seurat object (gitignored; loaded via LoadDataset).
## ============================================================================

source('03-scripts/R/seq_functions.R') 

if (!exists('nuclei')) { nuclei <- LoadDataset("Dec2024") }

subclass_colors <- LoadAllenColors('subclass')

SAVE_PLOTS <- FALSE
save_dir   <- "05-results/archive/neurogenesis_stage_axis"

ASSAY <- "SCT"   # per decision; set to "SCT" to mirror codebase default


# 1. Subclasses along the neurogenic axis (+ controls) -------------------------
# Ordered stem -> mature. Astro and Oligo are negative controls; CA1 is a
# non-DG mature-excitatory comparator. Matches the subclass set already used in
# the original SuppFigure7.R.
#   319 Astro-TE NN   : astroglia (shares RGL markers Gfap/Fabp7/Sox2 — control)
#   038 DG-PIR Ex IMN : immature-neuron subclass (expected to carry immature end)
#   037 DG Glut       : mature dentate granule cells
#   016 CA1-ProS Glut : mature non-DG excitatory comparator
#   327 Oligo NN      : oligodendrocyte negative control
neurogenic_subclasses <- c(
  "319 Astro-TE NN",
  "327 Oligo NN",
  "038 DG-PIR Ex IMN",
  "037 DG Glut"
)


# 2. A-priori gene panel (Nicola 2015 Fig 1A staging + standard extensions) -----
# `core` flags genes appearing directly in Nicola 2015 Fig 1A; extensions are
# widely-used stage markers cited to their own primary literature. Symbols are
# mouse MGI official symbols; common aliases noted.
stage_levels <- c("Type-1 RGL", "Type-2a IPC", "Type-2b/3 immature", "Mature GC")

neurogenesis_panel <- tibble::tribble(
  ~gene,      ~stage,                 ~core,  ~alias,         ~ref,
  # ---- Type-1 radial-glia-like / quiescent stem ----
  "Gfap",     "Type-1 RGL",           TRUE,   NA,             "Nicola 2015 Fig1A",
  "Id4",      "Type-1 RGL",           FALSE,  NA,             "Blomfield 2019; Zhang 2019",
  "Sox2",     "Type-1 RGL",           TRUE,   NA,             "Nicola 2015 Fig1A",
  "Hopx",     "Type-1 RGL",           FALSE,  NA,             "Berg 2019; Shin 2015",
  "Fabp7",    "Type-1 RGL",           TRUE,   "BLBP",         "Nicola 2015 (BLBP)",
  "Nes",      "Type-1 RGL",           TRUE,   "Nestin",       "Nicola 2015 Fig1A",
  # ---- Type-2a amplifying intermediate progenitor ----
  "Sox4",     "Type-2a IPC",          FALSE,  NA,             "Mu 2012",
  "Eomes",    "Type-2a IPC",          TRUE,   "Tbr2",         "Hodge 2008; Nicola 2015",
  "Mki67",    "Type-2a IPC",          TRUE,   "Ki67",         "Nicola 2015 (proliferation)",
  "Ascl1",    "Type-2a IPC",          FALSE,  "Mash1",        "Kim 2011; Andersen 2014",
  # ---- Type-2b / Type-3 neuronally-determined immature ----
  "Prox1",    "Type-2b/3 immature",   TRUE,   NA,             "Nicola 2015; Lavado 2010",
  "Dcx",      "Type-2b/3 immature",   TRUE,   "Doublecortin", "Nicola 2015 Fig1A",
  "Sox11",    "Type-2b/3 immature",   FALSE,  NA,             "Mu 2012; Haslinger 2009",
  "Igfbpl1",  "Type-2b/3 immature",   FALSE,  NA,             "Hagihara 2019",
  "Neurod1",  "Type-2b/3 immature",   TRUE,   "NeuroD1",      "Nicola 2015 Fig1A",
  "Calb2",    "Type-2b/3 immature",   FALSE,  "Calretinin",   "Brandt 2003 (CR immature)",
  # ---- Mature granule cell endpoint ----
  "Rbfox3",   "Mature GC",            TRUE,   "NeuN",         "Nicola 2015 (NeuN endpoint)",
  "Calb1",    "Mature GC",            FALSE,  "Calbindin",    "Brandt 2003 (mature GC)"
) |>
  mutate(
    stage      = factor(stage, levels = stage_levels),
    gene_order = row_number(),
    gene       = factor(gene, levels = gene)
  )


# 3. Subset to SE baseline across the neurogenic + control subclasses ----------
nuclei_ng <- nuclei |>
  subset((subclass_name %in% neurogenic_subclasses) & activity_condition == 'SE')

# fix subclass factor order to the stem -> mature axis
nuclei_ng$subclass_name <- factor(nuclei_ng$subclass_name, levels = neurogenic_subclasses)

# report what is actually present (sanity, not assumption)
print(table(nuclei_ng$subclass_name))

# ensure log-normalized data exists for the mean-expression metric on the RNA assay
DefaultAssay(nuclei_ng) <- ASSAY
if (ASSAY == "RNA") {
  nuclei_ng <- NormalizeData(nuclei_ng, assay = "RNA", verbose = FALSE)
}


# 4. Detection metrics (per gene x subclass) -----------------------------------
genes_all  <- as.character(neurogenesis_panel$gene)
counts_mat <- GetAssayData(nuclei_ng, assay = ASSAY, layer = "counts")
data_mat   <- GetAssayData(nuclei_ng, assay = ASSAY, layer = "data")
subv       <- nuclei_ng$subclass_name

genes_in   <- intersect(genes_all, rownames(counts_mat))
genes_out  <- setdiff(genes_all, rownames(counts_mat))
if (length(genes_out) > 0) {
  message("Panel genes ABSENT from ", ASSAY, " assay (drawn as not-in-assay): ",
          paste(genes_out, collapse = ", "))
}

detect_one <- function(g, sc) {
  cells <- which(subv == sc)
  n     <- length(cells)
  if (!g %in% rownames(counts_mat) || n == 0) {
    return(tibble(gene = g, subclass_name = sc, n_cells = n,
                  det_frac = NA_real_, mean_exp = NA_real_,
                  in_assay = g %in% rownames(counts_mat)))
  }
  tibble(
    gene          = g,
    subclass_name = sc,
    n_cells       = n,
    det_frac      = mean(counts_mat[g, cells] > 0),
    mean_exp      = mean(data_mat[g, cells]),
    in_assay      = TRUE
  )
}

detection <- tidyr::expand_grid(gene = genes_all,
                                subclass_name = levels(subv)) |>
  purrr::pmap_dfr(\(gene, subclass_name) detect_one(gene, subclass_name)) |>
  left_join(
    neurogenesis_panel |> mutate(gene = as.character(gene)) |>
      select(gene, stage, core, alias, ref, gene_order),
    by = "gene"
  ) |>
  arrange(gene_order) |>
  mutate(
    gene          = forcats::fct_reorder(gene, gene_order),
    subclass_name = factor(subclass_name, levels = neurogenic_subclasses)
  )

# write the auditable underlying table
if (SAVE_PLOTS) {
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(detection, "05-results/archive/neurogenesis_stage_axis/detection_table.csv")
}


# 5. Detection dotplot along the stage axis ------------------------------------
# x = subclass (stem -> mature, controls included)
# y = gene (ordered along Nicola axis, faceted by stage)
# size  = detection fraction; color = mean log-norm expression
# genes absent from the assay drawn as grey crosses
plot_df <- detection |>
  mutate(
    subclass_name = factor(subclass_name, levels = neurogenic_subclasses),
    absent        = !in_assay | is.na(det_frac)
  )

p_detect <- ggplot(plot_df, aes(x = subclass_name, y = gene)) +
  geom_point(data = ~ dplyr::filter(.x, !absent),
             aes(size = det_frac, color = mean_exp)) +
  scale_size_continuous(name = "Detection\nfraction",
                        limits = c(0, 1), range = c(0.3, 5)) +
  scale_color_viridis_c(name = "Mean log-norm\nexpression", option = "C") +
  facet_grid(stage ~ ., scales = "free_y", space = "free_y", switch = "y",
    labeller  = labeller(stage = c(
        "Type-1 RGL"         = "Type-1 RGL",
        "Type-2a IPC"        = "Type-2a IPC",
        "Type-2b/3 immature" = "Type-2b/3 immature",
        "Mature GC"          = "Mature\nGC"
      ))) +
  scale_y_discrete(expand = expansion(add = 0.75)) +
  labs(
    x = NULL, y = NULL,
  ) +
  theme_minimal(base_family = 'Helvetica') +
  theme(
    axis.text.x        = element_text(angle = 45, hjust = 1),
    axis.text.y        = element_text(face = "italic"),
    panel.grid.minor   = element_blank(),
    strip.text.y.left  = element_text(face = "bold", angle = 90),
    strip.placement    = "outside",
    legend.position    = "right"
  )

p_detect


# 6. Console summary — detection gradient by stage (review, not a claim) --------
summary_by_stage <- detection |>
  group_by(stage, subclass_name) |>
  summarise(mean_det_frac          = mean(det_frac, na.rm = TRUE),
            genes_detected_gt10pct = sum(det_frac > 0.10, na.rm = TRUE),
            n_genes                = dplyr::n(),
            .groups = "drop")
print(as.data.frame(summary_by_stage))


# 7. Save ----------------------------------------------------------------------
if (SAVE_PLOTS) {
  ggsave(glue("{save_dir}/detection_neurogenesis_stage_axis.png"),
         plot = p_detect, width = 4, height = 8, dpi = 300, bg = 'white')
  svgsave(plot = p_detect,
          filename = "detection_neurogenesis_stage_axis.svg",
          path = save_dir, width = 4, height = 8)
}
## ----
