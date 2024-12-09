# load libs and data ----
print("Loading libraries and data...")
library(Seurat)
library(tidyverse)
library(readxl)
library(patchwork)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)

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

# loop thru CSVs and load in DEGs
df_all_DEGs <- data.frame()  # init empty df
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

df_gene_classifications <- df_all_DEGs |> 
  mutate(is_ERG = ifelse(contrast == "d30m_vs_dSE", TRUE, FALSE),
         is_LRG = ifelse(contrast == "d6h_vs_dSE", TRUE, FALSE)) |> 
  group_by(gene) |> 
  summarize(
    is_ERG = any(is_ERG),
    is_LRG = any(is_LRG),
  ) |> 
  mutate(classification = case_when(
    is_ERG & is_LRG ~ "both",
    is_ERG ~ "ERG",
    is_LRG ~ "LRG",
    TRUE ~ NA_character_
  )) |> print()
# for each subclass, get expression for each condition

df_distinct_DEGs <- df_all_DEGs |> 
  mutate(condition = case_when(
    contrast == "d30m_vs_dSE" ~ "d30m",
    contrast == "d6h_vs_dSE" ~ "d6h",
    contrast == "KA_vs_dSE" ~ "KA",
    TRUE ~ NA_character_
  )) |>
  # filter out any contrasts that include KA
  filter(contrast != "KA_vs_dSE") |>
  group_by(gene, subclass, contrast) |> 
  summarize(avg_log2FC = avg_log2FC, .groups = 'drop') |> 
  left_join(df_gene_classifications, by = 'gene') |> 
  # Create new columns for ERG and LRG classification based on the contrasts
  mutate(subclass_by_contrast = paste(subclass, contrast, sep = " x ")) |> 
  mutate(n_subclasses = n_distinct(subclass), .by = 'gene') |> 
  distinct(gene, .keep_all = TRUE) |> 
  select(-subclass) # remove subclass for future joins

# normalize data ----
print("Normalizing expression data...")

# calculate log2FC expression compared to SE
nuclei_subclass <- subset(nuclei, subclass_name %in% subclass_list)

df_expression <- AverageExpression(
    nuclei_subclass,
    features = df_distinct_DEGs$gene,
    group.by = c("subclass_name", "condition"),
    assays = 'SCT',
    layer = 'data'
  )$SCT |> 
  as.data.frame() |> 
  rownames_to_column(var = 'gene') |> 
  pivot_longer(cols = -gene, names_to = 'subclass_by_contrast', values_to = 'avg_expression_ln_scale') |> 
  separate(subclass_by_contrast, into = c('subclass', 'condition'), sep = '_') |> 
  mutate(subclass = str_sub(subclass, 2, -1)) |> 
  left_join(df_distinct_DEGs, by = 'gene') |> 
  group_by(gene, subclass) |> 
  mutate(log2fc_from_SE = 
    (avg_expression_ln_scale - avg_expression_ln_scale[condition == "dSE"]) / log(2),
    .after = avg_expression_ln_scale
  ) |> 
  mutate(subclass_by_condition = paste(subclass, condition, sep = ' x ')) |> 
  mutate(gene = factor(gene, levels = df_distinct_DEGs$gene)) |> 
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
  arrange(desc(n_subclasses), subclass)
df_gene_levels_LRG <- df_expression |> 
  filter(classification == 'LRG') |> 
  filter(condition == 'd6h') |>
  slice_max(order_by = log2fc_from_SE, by = gene) |> 
  arrange(desc(n_subclasses), subclass)
df_gene_levels_both <- df_expression |> 
  filter(classification == 'both') |> 
  filter(condition == 'd30m' | condition == 'd6h') |>
  slice_max(order_by = log2fc_from_SE, by = gene) |> 
  arrange(desc(n_subclasses), subclass)

df_gene_levels <- rbind(df_gene_levels_ERG, df_gene_levels_LRG, df_gene_levels_both) |> 
  arrange(classification, desc(n_subclasses), subclass, desc(log2fc_from_SE))

# determine whether genes are TFs
mouse_TFs <- read_excel("02-data/published_data/Mus_musculus_TF.xlsx") |> 
  select(Symbol) |> 
  pull()

df_gene_levels <- df_gene_levels |> 
  mutate(TF = ifelse(gene %in% mouse_TFs, "1", "0"))

# combine classification and TF for anno purposes
df_gene_levels$classification_TF <- paste(df_gene_levels$classification, df_gene_levels$TF, sep = "_")

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


# sort genes and subclasses ----
# Create a matrix of log2 fold changes for the heatmap
expression_matrix <- df_expression |>
  filter(condition != 'dSE') |> 
  select(gene, subclass_by_condition, log2fc_from_SE) |>
  pivot_wider(names_from = gene, values_from = log2fc_from_SE) |> 
  column_to_rownames(var = "subclass_by_condition") |> 
  # select(-subclass) |> 
  as.matrix()

# re-order expression matrix
subclass_colors <- LoadAllenColors("subclass")
condition_colors <- LoadConditionColors("May2024")

# Ensure df_subclass_annotation has the same order as the rows in expression_matrix
df_subclass_annotation <- df_expression |>
  filter(condition != 'dSE') |> 
  group_by(subclass_by_condition, condition, subclass) |>
  distinct(subclass_by_condition) |> 
  arrange(condition, subclass)

# order expression matrix rows and cols
expression_matrix <- expression_matrix[match(df_subclass_annotation$subclass_by_condition, rownames(expression_matrix)), ]
expression_matrix <- expression_matrix[,match(df_gene_levels$gene, colnames(expression_matrix))]

# turn gene list (cols of exp matrix) into binary TF/non-TF
df_tf <- colnames(expression_matrix) |>
  as_tibble() |> 
  mutate(gene = colnames(expression_matrix)) |> 
  select(gene) |> 
  mutate(TF = ifelse(gene %in% mouse_TFs, "1", "0"))


# make annotation objects ----
row_anno <- rowAnnotation(
  condition = df_subclass_annotation$condition,
  subclass = df_subclass_annotation$subclass,
  col = list(
    condition = condition_colors,
    subclass = subclass_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

col_anno_top <- HeatmapAnnotation(
  classification = df_gene_levels$classification_TF,
  # subclass = df_gene_levels$subclass,
  col = list(
    classification = c(
      "ERG_0" = "#BB4430", "LRG_0" = "#F2B880", "both_0" = "#82A6B1", 
      "ERG_1" = "black", "LRG_1" = "black", "both_1" = "black"
      ),
    subclass = subclass_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

ieg_symbols <- LoadGeneList("IEG")
col_anno_bottom <- HeatmapAnnotation(
  iegs = anno_mark(
    at = which(df_gene_levels$gene %in% ieg_symbols), 
    labels = ieg_symbols, 
    side = "bottom"
    ),
  show_annotation_name = FALSE
)


# plot ----
p <- Heatmap(
  expression_matrix, # exclude the first column (subclass_by_condition)
  name = "log2FC",
  col = colorRamp2(c(-5,0,5), hcl_palette = "blue-red2"),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = row_anno,
  top_annotation = col_anno_top,
  bottom_annotation = col_anno_bottom,
  row_split = df_subclass_annotation$condition,
  column_split = factor(df_gene_levels$classification, levels = c("ERG", "LRG", "both")),
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
)
draw(p)
htShiny(p, action = 'hover', output_ui_float = T)


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