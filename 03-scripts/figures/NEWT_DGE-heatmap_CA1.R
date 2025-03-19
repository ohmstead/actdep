# This code produces a matrix heatmap of gene expression in 30m EE, 6h EE, and seizure compared to SE. The output is a 3-column matrix, where the columns are cells whose colormap corresponds to log2(FC) from SE expression. The columns are 30m EE, 6h EE, and Seizure. Each row is a DEG.

# - The rows of the matrix are sorted first according to ASCERTAINMENT, i.e. whether the gene was differentially expressed in 30m EE, 6h EE, or both.
# - Within each ascertainment group, the rows are sorted according to the PEAK EXPRESSION in either 30m EE or 6h EE.

# ---- load libs & data ----
print("Loading libraries and data...")

library(Seurat)
library(ggplot2)
library(readr)

source("03-scripts/R/seq_functions.R")
nuclei <- LoadDataset("May2024", "combined")

# ---- get DEG list ----
print("Getting DEG list...")
# read in genes from csv files containing previously calculated DEG results
dSE_v_d30m <- read_csv("04-analysis/DEGs/condition/016_CA1-ProS_Glut_d30m_vs_dSE.csv")
dSE_v_d6h <-  read_csv("04-analysis/DEGs/condition/016_CA1-ProS_Glut_d6h_vs_dSE.csv")
dSE_v_KA <-   read_csv("04-analysis/DEGs/condition/016_CA1-ProS_Glut_KA_vs_dSE.csv")

# get lists from SE vs 30m and SE vs 6h
degs <- unique(c(dSE_v_d30m$gene, dSE_v_d6h$gene))
genes_30m   <- degs[ (degs %in% dSE_v_d30m$gene) & !(degs %in% dSE_v_d6h$gene) ]
genes_6h    <- degs[ (degs %in% dSE_v_d6h$gene)  & !(degs %in% dSE_v_d30m$gene) ]
genes_both  <- degs[ (degs %in% dSE_v_d30m$gene) & (degs %in% dSE_v_d6h$gene) ]
degs <- c(genes_30m, genes_6h, genes_both)

# finally, remove any genes that are sexually dimorphically expressed between M/F
yao_sex_genes <- read_csv("04-analysis/DEGs/Yao2023_sex_DEGs/016_CA1-ProS_Glut_DEGs.csv")
degs <- degs[ !(degs %in% yao_sex_genes$gene) ]


# ---- norm expression data ----
print("Normalizing expression data...")
# get the average expression for each DEG grouped by condition
deg_matrix <- nuclei |> 
  subset(subclass_name == '016 CA1-ProS Glut') |> 
  AverageExpression(assays = 'SCT', layer = 'counts', features = degs, group.by = 'condition') |> 
  as.data.frame()

# get log2FC expression values (FC compared to SE)
deg_matrix_normed <- deg_matrix
deg_matrix_normed[, 2:4] <- deg_matrix_normed[, 2:4] + 1 / (deg_matrix_normed[, 1] + 1) # no zero division!
deg_matrix_normed[, 2:4] <- log2(deg_matrix_normed[, 2:4] + 1)
deg_matrix_normed[, 1] <- 1

colnames(deg_matrix_normed) <- c('SE', '30 min EE', '6 hr EE', 'Seizure')

# ---- sort gene list ----
print("Sorting gene list...")
# Convert the matrix to a data frame
df_degs <- as.data.frame(deg_matrix_normed)

# Add a column for gene names
df_degs$gene <- rownames(df_degs)

# add column representing whether df_degs$gene is in genes_30m, genes_6h, or genes_both
df_degs$ascertainment <- ifelse(
  df_degs$gene %in% genes_30m, '30 min EE', 
  ifelse(
    df_degs$gene %in% genes_6h, '6 hr EE', 
    'both'  # if not in either, then it must be in genes_both
  )
)

# pivot longer for plotting
df_degs_long <- df_degs |> 
  pivot_longer(cols = c("SE", "30 min EE", "6 hr EE", "Seizure"), names_to = "condition", values_to = "log2FC_expression") |>
  group_by(gene) |>
  relocate(ascertainment, .after = last_col())

# determine the peak expression condition for each gene
condition_info <- df_degs_long |> 
  filter(condition %in% c("30 min EE", "6 hr EE")) |> 
  group_by(gene, ascertainment) |> 
  summarize(peak_condition = condition[which.max(log2FC_expression)],
            peak_expression = max(log2FC_expression),
            .groups = 'drop')

# sort genes by peak condition and then by peak expression level
sorted_genes <- condition_info |> 
  mutate(peak_condition = factor(peak_condition, levels = c("6 hr EE", "30 min EE"))) |>
  mutate(ascertainment = factor(ascertainment, levels = c("both", "6 hr EE", "30 min EE"))) |>
  ungroup() |>
  group_by(peak_condition) |>
  arrange(.by_group = TRUE, desc(peak_expression)) |> 
  pull(gene)

# Update the gene factor levels based on the sorted order
df_degs_long$gene <- factor(df_degs_long$gene, levels = sorted_genes)


# ---- plot ----
print("Plotting heatmap...")

p <- ggplot(df_degs_long, aes(x = condition, y = gene, fill = log2FC_expression)) +
  geom_tile() +
  scale_fill_gradient2(low = 'blue', mid = 'white', high = 'red', limits = c(-3,3), midpoint = 0, oob = scales::squish) +
  theme_minimal() +
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        # axis.text.y = element_text(angle = 0, hjust = 1),
        axis.title = element_blank(),
        panel.grid = element_blank(),
        legend.position = 'none'
  ) +
  scale_x_discrete(limits = c("30 min EE", "6 hr EE", "Seizure"), expand = c(0,0)) +
  scale_y_discrete(position = 'right', guide = guide_axis(check.overlap = TRUE))

# ---- save plot ----
if (SAVE_PLOTS) {
  print("Saving plot...")
  print(p)
  save_path <- "05-results/NEWT_ZUPP_QC"
  ggsave(path = save_path, filename = 'DGE_heatmap_CA1.png', width = 2.843, height = 5.641, dpi = 900, units = 'in')
} else {
  print("Plotting without saving...")
  print(p)
}


print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))