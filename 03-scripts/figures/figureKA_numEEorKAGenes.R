library(ggplot2)
library(ggvenn)
library(readr)
library(dplyr)
library(glue)
library(patchwork)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset('Dec2024')
activity_colors <- LoadActivityColors()
subclass_colors <- LoadSubclassColors()


# subset object and get genes ----------------------------------------
# A helper function that safely reads a CSV file. If the file doesn't exist,
# returns an empty tibble with the expected columns.
safe_read_csv <- function(filepath) {
  if (file.exists(filepath)) {
    read_csv(filepath, show_col_types = FALSE)
  } else {
    warning(glue("File not found: {filepath}. Using an empty data frame."))
    tibble(gene = character(), classification = character())
  }
}

# Example: load the subclasses to iterate over.
subclasses_raw <- LoadSubclassesToUse(nuclei)
subclasses <- subclasses_raw[-c(2,8,10:14)] |> print() # dont load small subclasses

df_all <- tibble(df_all <- tibble(
  subclass    = character(),
  time        = character(),
  direction   = character(),
  EE_specific = integer(),
  KA_specific = integer(),
  shared      = integer()
))


# loop thru subclasses ------------------------------------------
for (subclass in subclasses) {
  subclass_str <- ShrinkSubclassName(subclass)
  
  file_deg_30m_EE <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__EE30m_vs_SE.csv")
  file_deg_30m_KA <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__KA30m_vs_SE.csv")
  file_deg_6h_EE  <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__EE6h_vs_SE.csv")
  file_deg_6h_KA  <- glue("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/{subclass_str}__KA6h_vs_SE.csv")
  
  deg_30m_EE <- safe_read_csv(file_deg_30m_EE)
  deg_30m_KA <- safe_read_csv(file_deg_30m_KA)
  deg_6h_EE  <- safe_read_csv(file_deg_6h_EE)
  deg_6h_KA  <- safe_read_csv(file_deg_6h_KA)
  
  safe_filter <- function(df, direction) {
    if ("classification" %in% names(df)) {
      df |> filter(classification == direction) |> pull(gene)
    } else {
      character()
    }
  }
  
  # Create a tibble summarizing counts for 30m and 6h for up and down regulated genes.
  df_current_subclass <- tibble(
    subclass = subclass,
    time = c('30m', '30m', '6h', '6h'),
    direction = c('up', 'down', 'up', 'down'),
    EE_specific = c(
      length(setdiff(safe_filter(deg_30m_EE, 'upregulated'),
                     safe_filter(deg_30m_KA, 'upregulated'))),
      length(setdiff(safe_filter(deg_30m_EE, 'downregulated'),
                     safe_filter(deg_30m_KA, 'downregulated'))),
      length(setdiff(safe_filter(deg_6h_EE, 'upregulated'),
                     safe_filter(deg_6h_KA, 'upregulated'))),
      length(setdiff(safe_filter(deg_6h_EE, 'downregulated'),
                     safe_filter(deg_6h_KA, 'downregulated')))
    ),
    KA_specific = c(
      length(setdiff(safe_filter(deg_30m_KA, 'upregulated'),
                     safe_filter(deg_30m_EE, 'upregulated'))),
      length(setdiff(safe_filter(deg_30m_KA, 'downregulated'),
                     safe_filter(deg_30m_EE, 'downregulated'))),
      length(setdiff(safe_filter(deg_6h_KA, 'upregulated'),
                     safe_filter(deg_6h_EE, 'upregulated'))),
      length(setdiff(safe_filter(deg_6h_KA, 'downregulated'),
                     safe_filter(deg_6h_EE, 'downregulated')))
    ),
    shared = c(
      length(intersect(safe_filter(deg_30m_EE, 'upregulated'),
                       safe_filter(deg_30m_KA, 'upregulated'))),
      length(intersect(safe_filter(deg_30m_EE, 'downregulated'),
                       safe_filter(deg_30m_KA, 'downregulated'))),
      length(intersect(safe_filter(deg_6h_EE, 'upregulated'),
                       safe_filter(deg_6h_KA, 'upregulated'))),
      length(intersect(safe_filter(deg_6h_EE, 'downregulated'),
                       safe_filter(deg_6h_KA, 'downregulated')))
    )
  )
  
  # if df_current_subclass is empty, skip it
  if (nrow(df_current_subclass) == 0) {
    next
  }
  
  print(df_current_subclass)
  df_all <- bind_rows(df_all, df_current_subclass)
}


df_all_processed <- df_all |>
  pivot_longer(
    cols = c(EE_specific, shared, KA_specific),
    names_to = "group",
    values_to = "count"
  ) |>
  mutate(
    # for downregulated genes, use negative counts
    count = ifelse(direction == "down", -count, count),
    time = factor(time, levels = c("30m", "6h")),
    subclass = factor(subclass, levels = (subclasses)),
    group = factor(group, levels = c("KA_specific", "shared", "EE_specific"))
  ) |> 
  print()

df_30m <- df_all_processed |> filter(time == "30m")
df_6h  <- df_all_processed |> filter(time == "6h")

gene_category_colors <- c(EE_specific = 'gray20', shared = 'gray50', KA_specific = 'gray80')


# plot num DEGs ------------------------------------------------------
# 30m plot
ggplot(df_30m, aes(y = count, x = subclass, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = 'black') +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = gene_category_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_blank(),
  )
si(1000, 500)

# 6h plot
ggplot(df_6h, aes(y = count, x = subclass, fill = group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), color = 'black') +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = gene_category_colors) +
  theme_bw() +
  theme(legend.position = 'none',
        legend.title = element_blank(),
        legend.text = element_text(size = 12),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 20),
        axis.text.x = element_blank(),
  )
si(1000,500)