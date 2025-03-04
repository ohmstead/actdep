print("Loading libraries and data...")
library(Seurat)
library(tidyverse)
library(readxl)
library(glue)
library(DESeq2)

source('03-scripts/R/seq_functions.R')

# load and subset object
nuclei <- readRDS('04-analysis/Seurats/Dec2024/seurat.Rds')

print('Subsetting Seurat object...')
subclass_list <- LoadSubclassesToUse(nuclei, ascertainment = 'custom')
nuclei_subclass <- subset(nuclei, subclass_name %in% subclass_list)
Idents(nuclei_subclass) <- nuclei_subclass$subclass_name


# read DEGs ----
print("Getting DEGs for all subclasses and contrasts...")
contrast_list <- c("EE30m_vs_SE", "EE6h_vs_SE")
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
  } else if (!(subclass %in% subclass_list)) {
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


# initial classification ---- 
print("Filtering DEGs")
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
  filter(abs(log2FoldChange) >= 0.585) |> 
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
  distinct(gene, .keep_all = TRUE) |> 
  mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) |> 
  select(-subclass) # remove subclass for future joins


# normalize data ----
print("Normalizing expression data...")
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
#                fc.name = 'log2FoldChange',
#                base = 2) |>
#     rownames_to_column('gene') |>
#     mutate(subclass = subclass,
#            contrast = comparison_name) |>
#     select(gene, log2FoldChange, subclass, contrast) # Keep necessary columns
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

# calculate log2FoldChange for each gene x subclass combo
df_expression.DESeq <- df_norm |>
  filter(!str_detect(activity_condition, 'KA')) |>
  group_by(gene, subclass) |>
  mutate(log2FC_calculated = avg_expression_log2 - avg_expression_log2[activity_condition == "SE"]) |>
  left_join(df_distinct_DEGs, by = 'gene') |>
  group_by(gene, subclass) |>
  mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
  mutate(activity_condition = factor(activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))) |>
  mutate(subclass = factor(subclass, levels = subclass_list)) |>
  relocate(gene, subclass, activity_condition, log2FC_calculated) |>
  ungroup()


# for genes in "both" category, force ERG/LRG ----
df_expression <- df_expression.DESeq
expression_range <- df_expression |>
  group_by(gene, subclass) |>
  summarise(
    max_expr_condition = activity_condition[which.max(log2FC_calculated)],
    min_expr_condition = activity_condition[which.min(log2FC_calculated)]
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


# save classifications ----
df_expression |> 
  separate(subclass_by_contrast, into = c("subclass", "contrast"), sep = " x ") |> 
  relocate(subclass, .after = gene) |> 
  dplyr::rename(log2FoldChange.pseudobulk = log2FoldChange) |>
  distinct(gene, classification, log2FoldChange.pseudobulk, subclass, contrast, is_ERG, is_LRG, gene) |> 
  write_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/0_DEG_classifications.csv") |> 
  print()
