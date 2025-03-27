library(ggplot2)
library(ggvenn)
library(readr)
library(dplyr)
library(glue)
library(patchwork)
library(ggpattern)

source('03-scripts/R/seq_functions.R')

# nuclei <- LoadDataset('Dec2024')
activity_colors <- LoadActivityColors()
subclass_colors <- LoadAllenColors()
save_dir <- "05-results/CRAB/raw_R_plots"

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
subclasses_raw <- LoadSubclassesToUse()
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
    values_to = "count_raw"
  ) |>
  # log10 transform counts
  mutate(count_log = log10(count_raw)) |> 
  # set log10(0) == -Inf >>> log10(0) == 0
  mutate(count_log = ifelse(is.infinite(count_log), 0, count_log)) |>
  mutate(
    # for downregulated genes, use negative counts
    count_log = ifelse(direction == "down", -count_log, count_log),
    label_y = ifelse(direction == "up", count_log + 0.2, count_log - 0.2),
    time = factor(time, levels = c("30m", "6h")),
    subclass = factor(subclass, levels = (subclasses)),
    group = factor(group, levels = c("KA_specific", "shared", "EE_specific"))
  ) |> 
  print()

gene_category_colors <- c(EE_specific = 'gray20', shared = 'gray50', KA_specific = 'gray80')
direction_patterns <- c(up = 'none', down = 'stripe')


# plot num DEGs ------------------------------------------------------
# 30m plot
p30m <- df_all_processed |> 
  filter(time == "30m") |> 
ggplot() +
  aes(y = count_log, x = subclass, fill = group) +
  geom_bar(stat = "identity", , position = position_dodge(width = 0.9), color = 'black') +
  geom_text(aes(label = count_raw, y = label_y), position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = gene_category_colors) +
  labs(title = 'DEG number by stimuli at 30m', y = 'log10(DEG count)') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    plot.title= element_text(size = 25)
  )
print(p30m)

# 6h plot
p6h <- df_all_processed |> 
  filter(time == "6h") |> 
  ggplot() +
  aes(y = count_log, x = subclass, fill = group) +
  geom_bar(stat = "identity", , position = position_dodge(width = 0.9), color = 'black') +
  geom_text(aes(label = count_raw, y = label_y), position = position_dodge(width = 0.9)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = gene_category_colors) +
  labs(title = 'DEG number by stimuli at 6h', y = 'log10(DEG count)') +
  theme_bw() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 15, angle = 30, hjust = 1),
    plot.title= element_text(size = 25)
  )
print(p6h)


# save --------------------------------------------
if (SAVE_PLOTS) {
  h = 8
  w = 12
  # png
  ggsave(p30m,
         filename = 'numDEG_30m.png',
         path = save_dir,
         height = h, width = w)
  ggsave(p6h,
         filename = 'numDEG_6h.png',
         path = save_dir,
         height = h, width = w)
  
  # svg
  svgsave(p30m + LoadBarebonesTheme(ticks = 'y'),
          filename = 'numDEG_30m.svg',
          savedir = save_dir,
          h = h, w = w)
  svgsave(p6h + LoadBarebonesTheme(ticks = 'y'),
          filename = 'numDEG_6h.svg',
          savedir = save_dir,
          h = h, w = w)
}
