# ---- load libs and data ----
print("Loading libraries and data...")
library(Seurat)
library(tidyverse)
library(readxl)
library(patchwork)

source('03-scripts/R/seq_functions.R')
nuclei <- LoadDataset("May2024", "combined")
condition_colors <- LoadConditionColors("May2024")


# ---- get DEGs for all subclasses and contrasts ----
print("Getting DEGs for all subclasses and contrasts...")

# get list of subclasses to use
subclass_list <- LoadSubclassesToUse(nuclei)

# get list of contrasts to use
contrast_list <- c("d30m_vs_dSE", "d6h_vs_dSE", "KA_vs_dSE")

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


# ---- get ordered, non-repeating list of DEGs ----
print("Culling duplicate DEGs and sorting")

# for each subclass, get expression for each condition
df_all_DEGs <- df_all_DEGs |> 
  group_by(gene, subclass, contrast) |> 
  summarize(avg_log2FC = avg_log2FC, .groups = 'drop') |> 
  print()

df_all_DEGs_classified <- df_all_DEGs |> 
  # Create new columns for ERG and LRG classification based on the contrasts
  mutate(is_ERG = ifelse(contrast == "d30m_vs_dSE", TRUE, FALSE),
         is_LRG = ifelse(contrast == "d6h_vs_dSE", TRUE, FALSE)) |> 
  # Group by gene to classify as ERG, LRG, or both
  group_by(gene) %>%
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

df_all_DEGs


# ---- normalize expression data ----
print("Normalizing expression data...")

# calculate log2FC expression compared to SE
nuclei_subclass <- subset(nuclei, subclass_name %in% subclass_list)

expression_data <-  AverageExpression(
    nuclei_subclass,
    features = df_all_DEGs_classified$gene,
    group.by = c("subclass_name", "condition"),
    assays = 'SCT',
    layer = 'counts'
  )$SCT |> 
  as.data.frame() |> 
  # let's manipulate the data to make it easier to plot
  rownames_to_column(var = 'gene') |>
  pivot_longer(cols = -gene, names_to = 'subclass_by_contrast', values_to = 'avg_expression') |> 
  separate(subclass_by_contrast, into = c('subclass', 'condition'), sep = '_') |> 
  mutate(subclass = str_remove(subclass, '^.')) |> 
  group_by(gene, subclass) |> 
  mutate(log2fc_from_SE = log2((avg_expression+1) / (avg_expression[condition == 'dSE']+1))) |> 
  mutate(subclass_by_condition = paste(subclass, condition, sep = ' x ')) |> 
  mutate(gene = factor(gene, levels = df_all_DEGs_classified$gene))

# ---- re-level values for plotting ----
print("Re-leveling values for plotting...")
# reformat subclass_by_condition for plotting
lvls_subclass_by_condition <- expression_data |> 
  arrange(condition, subclass) |> 
  ungroup() |> 
  distinct(subclass_by_condition) |> 
  pull(subclass_by_condition)

expression_data <- expression_data |> 
  mutate(subclass_by_condition = factor(
    subclass_by_condition, 
    levels = lvls_subclass_by_condition)
    )


# ---- plot ----
print("Plotting heatmap...")

p <- expression_data |> 
  filter(condition != 'dSE') |>
ggplot() +
  aes(x = gene, y = subclass_by_condition) +
  geom_point(aes(color = log2fc_from_SE, size = log2fc_from_SE)) +
  scale_color_gradient2(low = 'blue', mid = 'white', high = 'red', midpoint = 0, limits = c(-3,3)) +
  scale_size_binned_area(limits = c(-5, 5), oob = scales::squish) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  )


# ---- save plot ----
if (SAVE_PLOTS) {
  print("Saving plot...")
  ggsave(path = '05-results/figure1/', filename = 'DGE-heatmap_all_subclasses.png', plot = p, width = 10, height = 10, dpi = 900)
} else {
  print("Plotting without saving...")
  print(p)
}

print(glue("Script {basename(sys.frame(1)$ofile)} complete!"))