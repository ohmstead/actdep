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
subclass_list <- LoadSubclassesToUse(nuclei, ascertainment = 'custom')
nuclei_subclass <- subset(nuclei, subclass_name %in% subclass_list)
Idents(nuclei_subclass) <- nuclei_subclass$subclass_name

# get list of contrasts to use
contrast_list <- c("EE30m_vs_SE", "EE6h_vs_SE")

# load all DEGs from every subclass x contrast
dir_deg <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
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
    arrange(desc(log2FoldChange)) |> 
    filter(classification != 'no_change')
    # filter(padj < 0.05)
  
  deg$subclass <- subclass
  deg$contrast <- contrast
  
  df_all_DEGs <- rbind(df_all_DEGs, deg)
}


# classify DEGs ----
print("Culling duplicate DEGs and sorting")

df_gene_classifications <- df_all_DEGs |> 
  mutate(is_ERG = contrast == "EE30m_vs_SE",
         is_LRG = contrast == 'EE6h_vs_SE') |> 
  group_by(gene) |>
  summarize(
    is_ERG = any(is_ERG, na.rm = T),
    is_LRG = any(is_LRG, na.rm = T),
  ) |> 
  mutate(classification = case_when(
    is_ERG & is_LRG ~ "both",
    is_ERG ~ "ERG",
    is_LRG ~ "LRG",
    TRUE ~ NA_character_
  )) |> print()

# for each subclass, get expression for each activity_condition
df_distinct_DEGs <- df_all_DEGs |> 
  filter(subclass %in% subclass_list) |> 
  mutate(activity_condition = case_when(
    contrast == "EE30m_vs_SE" ~ "EE30m",
    contrast == "EE6h_vs_SE" ~ "EE6h",
    contrast == "KA30m_vs_SE" ~ "KA30m",
    contrast == "KA6h_vs_SE" ~ "KA6h",
    TRUE ~ NA_character_
  )) |>
  filter(!str_detect(contrast, "KA")) |>
  group_by(gene, subclass, contrast) |> 
  summarize(log2FoldChange = log2FoldChange, .groups = 'drop') |> 
  left_join(df_gene_classifications, by = 'gene') |> 
  mutate(subclass_by_contrast = paste(subclass, contrast, sep = " x ")) |> 
  mutate(n_subclasses = n_distinct(subclass), .by = 'gene') |> 
  distinct(gene, .keep_all = TRUE) |> 
  mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) |> 
  select(-subclass) # remove subclass for future joins


# normalize data ----
print("Normalizing expression data...")

# calculate log2FC expression compared to SE

# Seurat log2FC ----
# comparisons <- list(
#   "EE30m_vs_SE" = c("EE30m", "SE"),
#   "EE6h_vs_SE" = c("EE6h", "SE")
# )
# 
# # Function to calculate log2FC for a given subclass and comparison
# calculate_fc <- function(subclass, comparison_name, ident1, ident2) {
#   nuclei_subclass |>
#     FoldChange(group.by = 'activity_condition',
#                ident.1 = ident1,
#                ident.2 = ident2,
#                subset.ident = subclass,
#                features = df_distinct_DEGs$gene,
#                fc.name = 'log2fc_from_SE',
#                base = 2) |>
#     rownames_to_column('gene') |>
#     mutate(subclass = subclass,
#            contrast = comparison_name) |>
#     select(gene, log2fc_from_SE, subclass, contrast) # Keep necessary columns
# }
# 
# # Iterate over all subclasses and comparisons, storing results in a single tibble
# df_fc <- expand_grid(subclass = subclass_list, comparison = names(comparisons)) |>
#   mutate(ident1 = map_chr(comparison, ~ comparisons[[.x]][1]),
#          ident2 = map_chr(comparison, ~ comparisons[[.x]][2])) |>
#   pmap_dfr(~ calculate_fc(..1, ..2, ..3, ..4)) |>
#   tibble()
# 
# df_expression.seurat <- df_fc |>
#   separate(contrast, into = c('group1', 'group2'), sep = '_vs_') |>
#   mutate(activity_condition = factor(group1, levels = c('EE30m', 'EE6h'))) |>
#   select(-c(group1, group2)) |>
#   mutate(subclass = factor(subclass, levels = subclass_list)) |>
#   mutate(gene = factor(gene, levels = df_distinct_DEGs$gene)) |>
#   left_join(df_distinct_DEGs, by = 'gene') |>
#   mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
#   print()
  

# DESeq log2FC ----
# Extract cell-level metadata and raw counts from your Seurat object:
pseudobulk_counts <- AggregateExpression(
    nuclei_subclass,
    assays = 'RNA',
    features = df_distinct_DEGs$gene,
    group.by = c("subclass_name", "activity_condition")
  )$RNA |>
  as.matrix()

coldata <- data.frame(
  group = colnames(pseudobulk_counts),
  stringsAsFactors = FALSE
)
coldata <- coldata |>
  separate(group, into = c("subclass", "activity_condition"), sep = "_", remove = F)
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
df_expression.DESeq <- df_norm |>
  filter(!str_detect(activity_condition, 'KA')) |>
  group_by(gene, subclass) |>
  mutate(log2fc_from_SE = avg_expression_log2 - avg_expression_log2[activity_condition == "SE"]) |>
  left_join(df_distinct_DEGs, by = 'gene') |>
  group_by(gene, subclass) |>
  mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
  mutate(activity_condition = factor(activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))) |>
  mutate(subclass = factor(subclass, levels = subclass_list)) |>
  relocate(gene, subclass, activity_condition, avg_expression_log2, log2fc_from_SE) |>
  ungroup()


# for genes in the "both" category, force ERG/LRG ----
df_expression <- df_expression.DESeq
expression_range <- df_expression |>
  group_by(gene, subclass) |>
  summarise(
    max_expr_condition = activity_condition[which.max(log2fc_from_SE)],
    min_expr_condition = activity_condition[which.min(log2fc_from_SE)]
  )

# re-classify as ERG or LRG
df_expression <- df_expression |>
left_join(expression_range, by = c("gene", "subclass")) |>
mutate(classification = if_else(
  classification == "both",
  if_else(
    direction == "up",
    if_else(max_expr_condition == "EE30m", "ERG", "LRG"), # if upregulated
    if_else(min_expr_condition == "EE30m", "ERG", "LRG"),   # if downregulated
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


# sort gene list ----
print("Sorting gene list...")

# 1. Begin with a tibble with columns gene, subclass, activity_condition, classification (ERG/LRG/both) and log2FC
# 2. Filter tibble to focus on the 30m activity_condition for ERGs, the 6h activity_condition for LRGs, or both activity_conditions for "Both."
# 3. Group rows by gene.
# 4. Determine the subclass that has the highest log2FC expression.
# 5. Arrange by subclass.
# 6. Pull the genes and use them as the new sorting order.

# df_gene_levels_ERG <- df_expression |>
#   filter(classification == 'ERG') |>
#   filter(activity_condition == 'EE30m') |>
#   separate(subclass_by_contrast, into = c("subclass.ascertainment", "contrast"), sep = " x ", remove = FALSE) |>
#   filter(subclass == subclass.ascertainment)
#   # slice_max(order_by = log2fc_from_SE, by = gene)
# 
# df_gene_levels_LRG <- df_expression |>
#   filter(classification == 'LRG') |>
#   filter(activity_condition == 'EE6h') |>
#   separate(subclass_by_contrast, into = c("subclass.ascertainment", "contrast"), sep = " x ", remove = FALSE) |>
#   filter(subclass == subclass.ascertainment)
#   # slice_max(order_by = log2fc_from_SE, by = gene)
# 
# df_gene_levels <- rbind(df_gene_levels_ERG, df_gene_levels_LRG) |>
#   mutate(subclass.ascertainment = factor(subclass.ascertainment, levels = subclass_list)) |>
#   arrange(classification, desc(n_subclasses), subclass.ascertainment, desc(log2fc_from_SE))

# get gene classifications from df_expression
df_gene_levels <- df_distinct_DEGs |>
  select(-classification) |> # remove old classification that includes 'both' category
  left_join(distinct(select(df_expression, gene, classification)))

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
  filter(activity_condition %in% c('EE30m', 'EE6h')) |> 
  select(gene, subclass_by_activity_condition, log2fc_from_SE) |>
  pivot_wider(names_from = gene, values_from = log2fc_from_SE) |> 
  column_to_rownames(var = "subclass_by_activity_condition") |> 
  as.matrix()
  # scale()  # z-score the columns

# re-order expression_matrix rows
df_subclass_annotationtation <- df_expression |>
  filter(activity_condition != 'SE') |> 
  group_by(subclass_by_activity_condition, activity_condition, subclass) |>
  distinct(subclass_by_activity_condition) |> 
  arrange(activity_condition, subclass)

expression_matrix <- expression_matrix[match(df_subclass_annotationtation$subclass_by_activity_condition, rownames(expression_matrix)), ]


# seriate gene columns ----
expression_matrix <- expression_matrix[,df_gene_levels$gene]

erg_matrix <- expression_matrix[str_detect(rownames(expression_matrix), 'EE30m'), df_gene_levels$classification == 'ERG']
lrg_matrix <- expression_matrix[str_detect(rownames(expression_matrix), 'EE6h'),  df_gene_levels$classification == 'LRG']

meths <- c(
  'HC_average'
  # 'OLO_average',
  # 'HC_average'
)
for (meth in meths) {
erg_order <- seriate(erg_matrix, method = 'Heatmap', seriation_method = meth) |> get_order(2)
lrg_order <- seriate(lrg_matrix, method = 'Heatmap', seriation_method = meth) |> get_order(2)
col_order <- c(erg_order, lrg_order+length(erg_order))

# reorganize cols of expression matrix according to col_order
expression_matrix <- expression_matrix[, names(col_order)]

# relevel/arrange df_gene_levels for annotation objects
df_gene_levels <- df_gene_levels |> 
  mutate(gene = factor(gene, levels = colnames(expression_matrix))) |> 
  arrange(gene)


# make annotation objects ----
left_annotation <- rowAnnotation(
  activity_condition = df_subclass_annotationtation$activity_condition,
  subclass = df_subclass_annotationtation$subclass,
  col = list(
    activity_condition = activity_colors,
    subclass = subclass_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

right_annotation <- rowAnnotation(
  subclass = df_subclass_annotationtation$subclass,
  col = list(
    activity_condition = activity_colors,
    subclass = subclass_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

top_annotation <- HeatmapAnnotation(
  classification = df_gene_levels$classification_TF,
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
bottom_annotation <- HeatmapAnnotation(
  iegs = anno_mark(
    at = which(df_gene_levels$gene %in% ieg_symbols),
    labels = intersect(colnames(expression_matrix), ieg_symbols), 
    side = "bottom"
  ),
  show_annotation_name = FALSE
)


# plot ----
p <- Heatmap(
  expression_matrix, # exclude the first column (subclass_by_activity_condition)
  name = meth,
  heatmap_legend_param = list(
    title_position = "topcenter",
    direction = "horizontal"
  ),
  col = circlize::colorRamp2(c(-1, 0, 1), hcl_palette = 'Blue-Red 2'),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = left_annotation,
  right_annotation = right_annotation,
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  row_split = df_subclass_annotationtation$activity_condition,
  column_split = factor(df_gene_levels$classification, levels = c("ERG", "LRG")),
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  show_row_names = FALSE,
  show_column_names = FALSE,
  use_raster = FALSE
)
draw(p, heatmap_legend_side = 'bottom')
}
htShiny(p, action = 'hover', output_ui_float = T)


# plot just excitatory neurons ----
expression_matrix_excitatory <- expression_matrix[str_detect(rownames(expression_matrix), 'Glut'), ]
df_subclass_annotationtation_excitatory <- df_subclass_annotationtation |> 
  filter(str_detect(subclass, 'Glut'))

expression_matrix_excitatory <- expression_matrix_excitatory[,df_gene_levels$gene]

erg_matrix <- expression_matrix_excitatory[str_detect(rownames(expression_matrix_excitatory), 'EE30m'), df_gene_levels$classification == 'ERG']
lrg_matrix <- expression_matrix_excitatory[str_detect(rownames(expression_matrix_excitatory), 'EE6h'),  df_gene_levels$classification == 'LRG']

seriate_method <- 'HC_average'
  
erg_order <- seriate(erg_matrix, method = 'Heatmap', seriation_method = meth) |> get_order(2)
lrg_order <- seriate(lrg_matrix, method = 'Heatmap', seriation_method = meth) |> get_order(2)
col_order <- c(erg_order, lrg_order+length(erg_order))

# reorganize cols of expression matrix according to col_order
expression_matrix_excitatory <- expression_matrix_excitatory[, names(col_order)]

# relevel/arrange df_gene_levels for annotation objects
df_gene_levels <- df_gene_levels |> 
  mutate(gene = factor(gene, levels = colnames(expression_matrix_excitatory))) |> 
  arrange(gene)


# make annotation objects ----
left_annotation <- rowAnnotation(
  activity_condition = df_subclass_annotationtation_excitatory$activity_condition,
  subclass = df_subclass_annotationtation_excitatory$subclass,
  col = list(
    activity_condition = activity_colors,
    subclass = subclass_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

right_annotation <- rowAnnotation(
  subclass = df_subclass_annotationtation_excitatory$subclass,
  col = list(
    activity_condition = activity_colors,
    subclass = subclass_colors
  ),
  show_annotation_name = FALSE,
  show_legend = FALSE
)

top_annotation <- HeatmapAnnotation(
  classification = df_gene_levels$classification_TF,
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
bottom_annotation <- HeatmapAnnotation(
  iegs = anno_mark(
    at = which(df_gene_levels$gene %in% ieg_symbols),
    labels = intersect(colnames(expression_matrix_excitatory), ieg_symbols), 
    side = "bottom"
  ),
  show_annotation_name = FALSE
)


# plot ----
p <- Heatmap(
  expression_matrix_excitatory, # exclude the first column (subclass_by_activity_condition)
  name = meth,
  heatmap_legend_param = list(
    title_position = "topcenter",
    direction = "horizontal"
  ),
  col = circlize::colorRamp2(c(-3, 0, 3), hcl_palette = 'Blue-Red 2'),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  left_annotation = left_annotation,
  right_annotation = right_annotation,
  top_annotation = top_annotation,
  bottom_annotation = bottom_annotation,
  row_split = df_subclass_annotationtation_excitatory$activity_condition,
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
  png("05-results/figure1/raw_R_plots/DGE-heatmap_all_subclasses_pseudobulk.png", width = 15, height = 5, units = "in", res = 900)
  draw(p, heatmap_legend_side = 'bottom')
  dev.off()
} else {
  print("Plotting without saving...")
  print(p)
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))