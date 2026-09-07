## ---- SuppFig2_Supp1
# Make venn diagrams for DEGs at 30m and 6h in both CA1 and DG
source('03-scripts/R/seq_functions.R')

library(ggvenn)

# load data --------------------------------
deg_sets <- list(
  CA1_EE_30m = read_csv("04-analysis/dge_method_comparison/voom/016_CA1-ProS_Glut__EE30m_vs_SE.csv"),
  CA1_EE_6h  = read_csv("04-analysis/dge_method_comparison/voom/016_CA1-ProS_Glut__EE6h_vs_SE.csv"),
  CA1_KA_30m = read_csv("04-analysis/dge_method_comparison/voom/016_CA1-ProS_Glut__KA30m_vs_SE.csv"),
  CA1_KA_6h  = read_csv("04-analysis/dge_method_comparison/voom/016_CA1-ProS_Glut__KA6h_vs_SE.csv"),
  DG_EE_30m = read_csv("04-analysis/dge_method_comparison/voom/037_DG_Glut__EE30m_vs_SE.csv"),
  DG_EE_6h  = read_csv("04-analysis/dge_method_comparison/voom/037_DG_Glut__EE6h_vs_SE.csv"),
  DG_KA_30m = read_csv("04-analysis/dge_method_comparison/voom/037_DG_Glut__KA30m_vs_SE.csv"),
  DG_KA_6h  = read_csv("04-analysis/dge_method_comparison/voom/037_DG_Glut__KA6h_vs_SE.csv")
)

# pull just the gene column for DEGs from each item in list
deg_sets <- lapply(deg_sets,
       function(df)
         df |>
           filter(classification != 'no_change') |>
           pull(gene)
       )


# make venn diagrams --------------------------------
# CA1 30m
CA1_30m <- ggvenn(
  deg_sets[c("CA1_EE_30m", "CA1_KA_30m")],
  fill_color = c("#D81B60", "#098154FF"),
  auto_scale = T,
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")

# DG 30m
DG_30m <- ggvenn(
  deg_sets[c("DG_EE_30m", "DG_KA_30m")],
  fill_color = c("#D81B60", "#098154FF"),
  auto_scale = T,
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")

# CA1 6h
CA1_6h <- ggvenn(
  deg_sets[c("CA1_EE_6h", "CA1_KA_6h")],
  fill_color = c("#FFC107", "#63a45e"),
  auto_scale = T,
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")

# DG 6h
DG_6h <- ggvenn(
  deg_sets[c("DG_EE_6h", "DG_KA_6h")],
  fill_color = c("#FFC107", "#63a45e"),
  auto_scale = T,
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")

p1 <- CA1_30m + CA1_6h + DG_30m + DG_6h + plot_layout(nrow=1)
print(p1)

# save plots
if (exists("SAVE_PLOTS") & SAVE_PLOTS==TRUE) {
save_dir <- "05-results/Figure2_SuppFig1/raw_R_plots"
# png
ggsave(plot = CA1_30m,
       filename = "venn__CA1_30m.png",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = CA1_6h,
       filename = "venn__CA1_6h.png",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_30m,
       filename = "venn__DG_30m.png",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_6h,
       filename = "venn__DG_6h.png",
       width = 4, height = 4,
       path = save_dir)
# svg
ggsave(plot = CA1_30m,
       filename = "venn_CA1_30m.svg",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = CA1_6h,
       filename = "venn_CA1_6h.svg",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_30m,
       filename = "venn_DG_30m.svg",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_6h,
       filename = "venn_DG_6h.svg",
       width = 4, height = 4,
       path = save_dir)
}
  ## ----