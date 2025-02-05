# load libs and data ----
print("Loading libraries and data...")
library(Seurat)
library(tidyverse)
library(readxl)
library(DESeq2)
library(patchwork)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)

source('03-scripts/R/seq_functions.R')
nuclei <- readRDS('04-analysis/Seurats/Dec2024/seurat.Rds')
activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# read in DEGs ----
print("Getting DEGs for all subclasses and contrasts...")

# get list of subclasses to use
subclass_list <- LoadSubclassesToUse(nuclei) |> 
  pull(subclass_name)

custom_subclass_list <- c(
  '016 CA1-ProS Glut',
  '025 CA2-FC-IG Glut',
  '017 CA3 Glut',
  '037 DG Glut',
  '023 SUB-ProS Glut',
  '031 CT SUB Glut',
  '033 NP SUB Glut',
  '046 Vip Gaba',
  '047 Sncg Gaba',
  '048 RHP-COA Ndnf Gaba',
  '049 Lamp5 Gaba',
  '050 Lamp5 Lhx6 Gaba',
  '051 Pvalb chandelier Gaba',
  '052 Pvalb Gaba',
  '053 Sst Gaba',
  '038 DG-PIR Ex IMN',
  '319 Astro-TE NN',
  '326 OPC NN',
  '327 Oligo NN',
  '334 Microglia NN'
)

# get list of contrasts to use
contrast_list <- c("EE30m_vs_SE", "EE6h_vs_SE")

# load all DEGs from every subclass x contrast
dir_deg <- "04-analysis/DEGs/Dec2024_activity_condition"
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
  subclass <- str_sub(subclass, 1, -2)
  # skip excluded subclasses or contrasts
  if (!(contrast %in% contrast_list)) {
    next
  }
  
  # add gene list to df
  deg <- read_csv(file) |> 
    arrange(desc(avg_log2FC))
  deg$subclass <- subclass
  deg$contrast <- contrast
  df_all_DEGs <- rbind(df_all_DEGs, deg)
}

df_all_DEGs <- df_all_DEGs |> 
  filter(p_val_adj < 0.05)


# classify DEGs ----
print("Culling duplicate DEGs and sorting")

df_gene_classifications <- df_all_DEGs |> 
  mutate(is_ERG = ifelse(contrast == "EE30m_vs_SE", TRUE, FALSE),
         is_LRG = ifelse(contrast == "EE6h_vs_SE", TRUE, FALSE)) |> 
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

# for each subclass, get expression for each activity_condition
df_distinct_DEGs <- df_all_DEGs |> 
  filter(subclass %in% custom_subclass_list) |> 
  mutate(activity_condition = case_when(
    contrast == "EE30m_vs_SE" ~ "EE30m",
    contrast == "EE6h_vs_SE" ~ "EE6h",
    contrast == "KA30m_vs_SE" ~ "KA30m",
    contrast == "KA6h_vs_SE" ~ "KA6h",
    TRUE ~ NA_character_
  )) |>
  filter(!str_detect(contrast, "KA")) |>
  group_by(gene, subclass, contrast) |> 
  summarize(avg_log2FC = avg_log2FC, .groups = 'drop') |> 
  left_join(df_gene_classifications, by = 'gene') |> 
  mutate(subclass_by_contrast = paste(subclass, contrast, sep = " x ")) |> 
  mutate(n_subclasses = n_distinct(subclass), .by = 'gene') |> 
  distinct(gene, .keep_all = TRUE) |> 
  mutate(direction = ifelse(avg_log2FC > 0, 'up', 'down')) |> 
  select(-subclass) # remove subclass for future joins


# normalize data ----
print("Normalizing expression data...")

# calculate log2FC expression compared to SE
nuclei_subclass <- subset(nuclei, subclass_name %in% custom_subclass_list)

# try DESeq2 approach instead
# Extract cell-level metadata and raw counts from your Seurat object:
pseudobulk_counts <- AggregateExpression(
    nuclei_subclass,
    assays = 'RNA',
    features = df_distinct_DEGs$gene,
    group.by = c("subclass_name", "activity_condition")
  )$RNA

coldata <- data.frame(
  group = colnames(pseudobulk_counts),
  stringsAsFactors = FALSE
)
coldata <- coldata |> 
  mutate(
    subclass = sapply(strsplit(group, "_"), `[`, 1),
    activity_condition = sapply(strsplit(group, "_"), `[`, 2)
  )
rownames(coldata) <- coldata$group

dds <- DESeqDataSetFromMatrix(
  countData = pseudobulk_counts,
  colData = coldata,
  design = ~ subclass + activity_condition
)

# Set "SE" as the reference level for activity_condition:
dds$activity_condition <- relevel(dds$activity_condition, ref = "SE")
dds <- DESeq(dds)

# log-transform normalized counts
norm_counts <- counts(dds, normalized = TRUE)
norm_log2 <- log2(norm_counts + 1)  # Adding a pseudocount to avoid log(0)
df_norm <- as.data.frame(norm_log2) |> 
  rownames_to_column(var = "gene") |> 
  pivot_longer(cols = -gene, names_to = "group", values_to = "avg_expression_log2")

# re-extract variables
df_norm <- df_norm |> 
  separate(group, into = c("subclass", "activity_condition"), sep = "_") |> 
  mutate(subclass = str_sub(subclass, 2, -1))

# calculate log2fc_from_SE for each gene x subclass combo
df_expression <- df_norm |> 
  group_by(gene, subclass) |> 
  mutate(log2fc_from_SE = avg_expression_log2 - avg_expression_log2[activity_condition == "SE"]) |> 
  left_join(df_distinct_DEGs, by = 'gene') |>
  group_by(gene, subclass) |> 
  mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |> 
  mutate(gene = factor(gene, levels = df_distinct_DEGs$gene)) |> 
  mutate(activity_condition = factor(activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))) |> 
  mutate(subclass = factor(subclass, levels = custom_subclass_list)) |> 
  relocate(gene, subclass, activity_condition, avg_expression_log2, log2fc_from_SE) |> 
  ungroup()


# for genes in the "both" category, refactor ----
expr_range <- df_expression |>
  filter(activity_condition == 'EE30m' | activity_condition == 'EE6h') |> 
  group_by(gene, subclass) |>
  summarise(
    max_expr_condition = activity_condition[which.max(avg_expression_log2)],
    min_expr_condition = activity_condition[which.min(avg_expression_log2)]
  )

# re-classify as ERG or LRG
df_expression <- df_expression |>
left_join(expr_range, by = c("gene", "subclass")) |>
mutate(classification = if_else(
  classification == "both",
  if_else(
    direction == "up",
    if_else(max_expr_condition == "EE30m", "ERG", "LRG"), # if upregulated
    if_else(min_expr_condition == "EE30m", "ERG", "LRG"   # if downregulated
      ),
    ),
  classification  # For genes not labeled "both", keep original classification.
  ),
)

# for any genes that may be ERG in one subclass and LRG in another, force to modal classification
df_expression <- df_expression |>
  group_by(gene) |>
  mutate(modal_classification = names(sort(table(classification), decreasing = TRUE))[1]) |>
  ungroup() |>
  mutate(classification = modal_classification) |>
  select(-modal_classification)


# take out seizure conditions----
df_expression <- df_expression |> 
  filter(!str_detect(activity_condition, 'KA'))

# sort gene list ----
print("Sorting gene list...")

# 1. Begin with a tibble with columns gene, subclass, activity_condition, classification (ERG/LRG/both) and log2FC
# 2. Filter tibble to focus on the 30m activity_condition for ERGs, the 6h activity_condition for LRGs, or both activity_conditions for "Both."
# 3. Group rows by gene.
# 4. Determine the subclass that has the highest log2FC expression.
# 5. Arrange by subclass.
# 6. Pull the genes and use them as the new sorting order.

df_gene_levels_ERG <- df_expression |> 
  filter(classification == 'ERG') |> 
  filter(activity_condition == 'EE30m') |>
  slice_max(order_by = log2fc_from_SE, by = gene) |> 
  arrange(desc(n_subclasses), subclass)
df_gene_levels_LRG <- df_expression |> 
  filter(classification == 'LRG') |> 
  filter(activity_condition == 'EE6h') |>
  slice_max(order_by = log2fc_from_SE, by = gene) |> 
  arrange(desc(n_subclasses), subclass)

# df_gene_levels <- rbind(df_gene_levels_ERG, df_gene_levels_LRG, df_gene_levels_both) |> 
df_gene_levels <- rbind(df_gene_levels_ERG, df_gene_levels_LRG) |> 
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
# reformat subclass_by_activity_condition for plotting
lvls_subclass_by_activity_condition <- df_expression |> 
  arrange(activity_condition, subclass) |> 
  ungroup() |> 
  distinct(subclass_by_activity_condition) |> 
  pull(subclass_by_activity_condition)

df_expression <- df_expression |> 
  mutate(subclass_by_activity_condition = factor(subclass_by_activity_condition, levels = lvls_subclass_by_activity_condition))


# sort genes and subclasses ----
# Create a matrix of log2 fold changes for the heatmap
expression_matrix <- df_expression |>
  filter(activity_condition != 'SE') |> 
  select(gene, subclass_by_activity_condition, log2fc_from_SE) |>
  pivot_wider(names_from = gene, values_from = log2fc_from_SE) |> 
  column_to_rownames(var = "subclass_by_activity_condition") |> 
  as.matrix() |> 
  scale()  # z-score the columns 

# re-order expression matrix
# Ensure df_subclass_annotation has the same order as the rows in expression_matrix
df_subclass_annotation <- df_expression |>
  filter(activity_condition != 'SE') |> 
  group_by(subclass_by_activity_condition, activity_condition, subclass) |>
  distinct(subclass_by_activity_condition) |> 
  arrange(activity_condition, subclass)

# order expression matrix rows and cols
expression_matrix <- expression_matrix[match(df_subclass_annotation$subclass_by_activity_condition, rownames(expression_matrix)), ]
expression_matrix <- expression_matrix[,match(df_gene_levels$gene, colnames(expression_matrix))]

# turn gene list (cols of exp matrix) into binary TF/non-TF
df_tf <- colnames(expression_matrix) |>
  as_tibble() |> 
  mutate(gene = colnames(expression_matrix)) |> 
  select(gene) |> 
  mutate(TF = ifelse(gene %in% mouse_TFs, "1", "0"))


# make annotation objects ----
left_anno <- rowAnnotation(
  activity_condition = df_subclass_annotation$activity_condition,
  subclass = df_subclass_annotation$subclass,
  col = list(
    activity_condition = activity_colors,
    subclass = subclass_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

right_anno <- rowAnnotation(
  subclass = df_subclass_annotation$subclass,
  col = list(
    activity_condition = activity_colors,
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
# re-arrange IEG symbols to match their order in the heatmap
iegs_ordered = df_gene_levels$gene[which(df_gene_levels$gene %in% ieg_symbols)]
col_anno_bottom <- HeatmapAnnotation(
  iegs = anno_mark(
    at = which(df_gene_levels$gene %in% ieg_symbols), 
    labels = iegs_ordered, 
    side = "bottom"
    ),
  show_annotation_name = FALSE
)


# plot ----
p <- Heatmap(
  expression_matrix, # exclude the first column (subclass_by_activity_condition)
  name = "z-scored log2FC from SE",
  heatmap_legend_param = list(
    title_position = "topcenter",
    at = c(-2, 0, 2),
    direction = "horizontal"
  ),
  col = circlize::colorRamp2(c(-2.5, 0, 2.5), hcl_palette = 'Blue-Red 2'),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = left_anno,
  right_annotation = right_anno,
  top_annotation = col_anno_top,
  bottom_annotation = col_anno_bottom,
  row_split = df_subclass_annotation$activity_condition,
  column_split = factor(df_gene_levels$classification, levels = c("ERG", "LRG")),
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = FALSE
)
draw(p, heatmap_legend_side = 'bottom')
htShiny(p, action = 'hover', output_ui_float = T)


# save plot ----
if (SAVE_PLOTS) {
  print("Saving plot...")
  png("05-results/figure1/raw_R_plots/DGE-heatmap_all_subclasses.png", width = 15, height = 5, units = "in", res = 900)
  draw(p, heatmap_legend_side = 'bottom')
  dev.off()
´} else {
  print("Plotting without saving...")
  print(p)
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))