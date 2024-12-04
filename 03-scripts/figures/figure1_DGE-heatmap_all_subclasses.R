# load libs and data ----
print("Loading libraries and data...")
library(Seurat)
library(tidyverse)
library(readxl)
library(patchwork)
library(ComplexHeatmap)

source('03-scripts/R/seq_functions.R')
nuclei <- LoadDataset("May2024", "combined")
condition_colors <- LoadConditionColors("May2024")


# read in DEGs ----
print("Getting DEGs for all subclasses and contrasts...")

# get list of subclasses to use
subclass_list <- LoadSubclassesToUse(nuclei)

# get list of contrasts to use
contrast_list <- c("d30m_vs_dSE", "d6h_vs_dSE")

# load all DEGs from every subclass x contrast
dir_deg <- "04-analysis/DEGs/condition/"
deg_files <- list.files(dir_deg, full.names = TRUE)

# make a df with columns for subclass and contrast
df_all_DEGs <- data.frame()  # init empty df

# loop thru CSVs and load in DEGs
for (file in deg_files) {
  fname <- basename(file)
  fname <- str_remove(fname, "\\.csv$")
  # Extract the contrast (last two underscore-separated parts in filename)
  contrast <- str_extract(fname, "[^_]+_[^_]+_[^_]+$")
  # Extract the subclass (everything before the contrast in filename)
  subclass <- str_remove(fname, paste0("_", contrast, "$"))
  subclass <- str_replace_all(subclass, "_", " ")
  # skip excluded subclasses or contrasts
  if (!(subclass %in% subclass_list)) {
    next
  }
  else if (!(contrast %in% contrast_list)) {
    next
  }
  
  # add gene list to df
  deg <- read_csv(file)
  deg$subclass <- subclass
  deg$contrast <- contrast
  df_all_DEGs <- rbind(df_all_DEGs, deg)
}


# classify DEGs ----
print("Culling duplicate DEGs and sorting")

# for each subclass, get expression for each condition
df_all_DEGs <- df_all_DEGs |> 
  # filter out any contrasts that include KA
  filter(contrast != "KA_vs_dSE") |>
  group_by(gene, subclass, contrast) |> 
  summarize(avg_log2FC = avg_log2FC, .groups = 'drop') |> 
  print()

df_all_DEGs_classified <- df_all_DEGs |> 
  # Create new columns for ERG and LRG classification based on the contrasts
  mutate(is_ERG = ifelse(contrast == "d30m_vs_dSE", TRUE, FALSE),
         is_LRG = ifelse(contrast == "d6h_vs_dSE", TRUE, FALSE)) |> 
  # Group by gene to classify as ERG, LRG, or both
  group_by(gene) |> 
  summarize(
    is_ERG = any(is_ERG),
    is_LRG = any(is_LRG),
    n_subclasses = n_distinct(subclass)
  ) |> 
  # Classify genes as ERG, LRG, or both
  mutate(classification = case_when(
    is_ERG & is_LRG ~ "both",
    is_ERG ~ "ERG",
    is_LRG ~ "LRG",
    TRUE ~ NA_character_
  )) |> 
  # Keep only distinct genes and remove rows without a classification
  filter(!is.na(classification)) |> 
  distinct(gene, .keep_all = TRUE) |> 
  # Arrange by classification and number of subclasses
  arrange(
    desc(classification == "ERG"), 
    desc(classification == "LRG"), 
    desc(classification == "both"), 
    desc(n_subclasses)
    )

df_all_DEGs <- df_all_DEGs |> 
  left_join(df_all_DEGs_classified, by = 'gene') |> 
  mutate(gene = factor(gene, levels = df_all_DEGs_classified$gene)) |>
  mutate(subclass_by_contrast = paste(subclass, contrast, sep = " x "))


# normalize data ----
print("Normalizing expression data...")

# calculate log2FC expression compared to SE
nuclei_subclass <- subset(nuclei, subclass_name %in% subclass_list)

df_expression <-  AverageExpression(
    nuclei_subclass,
    features = df_all_DEGs_classified$gene,
    group.by = c("subclass_name", "condition"),
    assays = 'SCT',
    layer = 'counts'
  )$SCT |> 
  as.data.frame() |> 
  rownames_to_column(var = 'gene') |>
  pivot_longer(cols = -gene, names_to = 'subclass_by_contrast', values_to = 'avg_expression') |> 
  separate(subclass_by_contrast, into = c('subclass', 'condition'), sep = '_') |> 
  mutate(subclass = str_remove(subclass, '^.')) |> 
  left_join(df_all_DEGs_classified, by = 'gene') |>
  group_by(gene, subclass) |> 
  mutate(log2fc_from_SE = log2((avg_expression+1) / (avg_expression[condition == 'dSE']+1))) |> 
  mutate(subclass_by_condition = paste(subclass, condition, sep = ' x ')) |> 
  mutate(gene = factor(gene, levels = df_all_DEGs_classified$gene)) |> 
  mutate(condition = factor(condition, levels = c('dSE', 'd30m', 'd6h', 'KA'))) |> 
  mutate(subclass = factor(subclass, levels = subclass_list)) |> 
  ungroup()


# sort gene list ----
print("Sorting gene list...")

# 1. Begin with a tibble with columns gene, subclass, condition, classification (ERG/LRG/both) and log2FC
# 2. Filter tibble to focus on the 30m condition for ERGs, the 6h condition for LRGs, or both conditions for "Both."
# 3. Group rows by gene.
# 4. Determine the subclass that has the highest log2FC expression.
# 5. Arrange by subclass.
# 6. Pull the genes and use them as the new sorting order.

df_gene_levels_ERG <- df_expression |> 
  filter(classification == 'ERG') |> 
  filter(condition == 'd30m') |>
  slice_max(order_by = log2fc_from_SE, by = gene) |> 
  arrange(desc(n_subclasses), subclass) |> print()
df_gene_levels_LRG <- df_expression |> 
  filter(classification == 'LRG') |> 
  filter(condition == 'd6h') |>
  slice_max(order_by = log2fc_from_SE, by = gene) |> 
  arrange(desc(n_subclasses), subclass) |> print()
df_gene_levels_both <- df_expression |> 
  filter(classification == 'both') |> 
  filter(condition == 'd30m' | condition == 'd6h') |>
  slice_max(order_by = log2fc_from_SE, by = gene) |> 
  arrange(desc(n_subclasses), subclass) |> print()
df_gene_levels <- rbind(df_gene_levels_ERG, df_gene_levels_LRG, df_gene_levels_both) |> 
  arrange(classification, desc(n_subclasses), subclass, desc(log2fc_from_SE)) |> print()

# re-level variables for plotting ----
print("Re-leveling variables for plotting...")
# reformat subclass_by_condition for plotting
lvls_subclass_by_condition <- df_expression |> 
  arrange(condition, subclass) |> 
  ungroup() |> 
  distinct(subclass_by_condition) |> 
  pull(subclass_by_condition)

df_expression <- df_expression |> 
  mutate(subclass_by_condition = factor(subclass_by_condition, levels = lvls_subclass_by_condition))


# sort genes and subclasses -------
# Create a matrix of log2 fold changes for the heatmap
expression_matrix <- df_expression |>
  filter(condition != 'dSE') |> 
  select(gene, subclass_by_condition, log2fc_from_SE) |>
  pivot_wider(names_from = gene, values_from = log2fc_from_SE) |> 
  column_to_rownames(var = "subclass_by_condition") |> 
  # select(-subclass) |> 
  as.matrix()

# make annotations ----
subclass_colors <- LoadAllenColors("subclass")
condition_colors <- LoadConditionColors("May2024")

# Ensure subclass_annotation has the same order as the rows in expression_matrix
subclass_annotation <- df_expression |>
  filter(condition != 'dSE') |> 
  group_by(subclass_by_condition, condition, subclass) |>
  distinct(subclass_by_condition) |> 
  arrange(condition, subclass)

# Define row annotation
row_anno <- rowAnnotation(
  condition = subclass_annotation$condition,
  subclass = subclass_annotation$subclass,
  col = list(
    condition = condition_colors,
    subclass = subclass_colors
  )
)

# Define col annotation
col_anno <- HeatmapAnnotation(
  classification = df_gene_levels$classification,
  # subclass = df_gene_levels$subclass,
  col = list(
    classification = c("ERG" = "#BB4430", "LRG" = "#F2B880", "both" = "#82A6B1"),
    subclass = subclass_colors
  )
)

# order expression matrix rows and cols
expression_matrix <- expression_matrix[match(subclass_annotation$subclass_by_condition, rownames(expression_matrix)), ]
expression_matrix <- expression_matrix[,match(df_gene_levels$gene, colnames(expression_matrix))]

# plot ----
p <- Heatmap(
  expression_matrix, # exclude the first column (subclass_by_condition)
  name = "log2FC",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = row_anno,
  top_annotation = col_anno,
  row_split = subclass_annotation$condition,
  column_split = factor(df_gene_levels$classification, levels = c("ERG", "LRG", "both")),
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  # row_names_gp = gpar(fontsize = 8),
  # column_names_gp = gpar(fontsize = 8, rot = 90)
)


# save plot ----
if (SAVE_PLOTS) {
  print("Saving plot...")
  # ggsave(path = '05-results/figure1/raw_R_plots/', filename = 'DGE-heatmap_all_subclasses.png', plot = p, width = 15, height = 5, dpi = 900)
  png("05-results/figure1/raw_R_plots/DGE-heatmap_all_subclasses.png", width = 15, height = 5, units = "in", res = 900)
  print(p)
  dev.off()
} else {
  print("Plotting without saving...")
  print(p)
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))