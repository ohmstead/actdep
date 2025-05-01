# Make venn diagrams for DEGs at 30m and 6h in both CA1 and DG
source('03-scripts/R/seq_functions.R')

library(ggvenn)

# load data --------------------------------
deg_sets <- list(
  CA1_EE_30m = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/016_CA1-ProS_Glut__EE30m_vs_SE.csv"),
  CA1_EE_6h  = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/016_CA1-ProS_Glut__EE6h_vs_SE.csv"),
  CA1_KA_30m = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/016_CA1-ProS_Glut__KA30m_vs_SE.csv"),
  CA1_KA_6h  = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/016_CA1-ProS_Glut__KA6h_vs_SE.csv"),
  DG_EE_30m = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/037_DG_Glut__EE30m_vs_SE.csv"),
  DG_EE_6h  = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/037_DG_Glut__EE6h_vs_SE.csv"),
  DG_KA_30m = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/037_DG_Glut__KA30m_vs_SE.csv"),
  DG_KA_6h  = read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/037_DG_Glut__KA6h_vs_SE.csv")
)

# pull just the gene column for DEGs from each item in list
deg_sets <- lapply(sets_degs, 
       function(df) 
         df |> 
         filter(padj < 0.05 & abs(log2FoldChange.raw) > 0.585) |> 
         pull(gene)
       )


# make venn diagrams --------------------------------
# CA1 30m
CA1_30m <- ggvenn(
  deg_sets[c("CA1_EE_30m", "CA1_KA_30m")],
  fill_color = c("#D81B60", "#098154FF"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")
CA1_30m

# CA1 6h
CA1_6h <- ggvenn(
  deg_sets[c("CA1_EE_6h", "CA1_KA_6h")],
  fill_color = c("#D81B60", "#098154FF"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")
CA1_6h

# DG 30m
DG_30m <- ggvenn(
  deg_sets[c("DG_EE_30m", "DG_KA_30m")],
  fill_color = c("#D81B60", "#098154FF"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")
DG_30m

# DG 6h
DG_6h <- ggvenn(
  deg_sets[c("DG_EE_6h", "DG_KA_6h")],
  fill_color = c("#D81B60", "#098154FF"),
  stroke_size = 0.5,
  set_name_size = 4,
  show_percentage = FALSE
) +
  theme(legend.position = "none")
DG_6h

# save plots
save_dir <- "05-results/CRAB/raw_R_plots"
# png
ggsave(plot = CA1_30m,
       filename = "05-results/CRAB/raw_R_plots/venn_CA1_30m.png",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = CA1_6h,
       filename = "05-results/CRAB/raw_R_plots/venn_CA1_6h.png",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_30m,
       filename = "05-results/CRAB/raw_R_plots/venn_DG_30m.png",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_6h,
       filename = "05-results/CRAB/raw_R_plots/venn_DG_6h.png",
       width = 4, height = 4,
       path = save_dir)
# svg
ggsave(plot = CA1_30m,
       filename = "05-results/CRAB/raw_R_plots/venn_CA1_30m.svg",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = CA1_6h,
       filename = "05-results/CRAB/raw_R_plots/venn_CA1_6h.svg",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_30m,
       filename = "05-results/CRAB/raw_R_plots/venn_DG_30m.svg",
       width = 4, height = 4,
       path = save_dir)
ggsave(plot = DG_6h,
       filename = "05-results/CRAB/raw_R_plots/venn_DG_6h.svg",
       width = 4, height = 4,
       path = save_dir)
