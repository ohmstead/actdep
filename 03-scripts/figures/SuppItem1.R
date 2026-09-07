# Supplementary Item 1: per-replicate metadata table (csv)
#
# Builds one row per biological replicate with the same columns shown in
# 05-results/SuppItem1/all_samples.png (activity_condition, age, sex,
# hemisphere, ZT, nuclei) plus one nuclei-count column for each of the 20
# subclasses returned by LoadSubclassesToUse().

source('03-scripts/R/seq_functions.R')

df_meta <- LoadMeta()
subclass_list <- LoadSubclassesToUse()

# seurat_meta.rds is missing age/hemisphere for the ZT12_30m replicates, so
# backfill those from the documentation csv the same way SuppTable1.R does.
df_metadata_csv <- read_csv(
  "01-documentation/sample_experimental_metadata_Dec2024.csv",
  show_col_types = FALSE
)

# per-replicate metadata + total nuclei ------------------------------------
df_sample <- df_meta |>
  as_tibble() |>
  select(sample, age, sex, hemisphere_seq, ZT, activity_condition) |>
  left_join(df_metadata_csv, by = 'sample', suffix = c('.seurat', '')) |>
  group_by(sample) |>
  summarize(
    age                = first(age),
    sex                = first(sex),
    hemisphere         = first(hemisphere_seq),
    ZT                 = first(ZT),
    activity_condition = first(activity_condition),
    nuclei             = n(),
    .groups = 'drop'
  ) |>
  arrange(activity_condition, sample) |>
  mutate(sample_num = seq_along(sample)) |>
  select(sample_num, sample, activity_condition, age, sex, hemisphere, ZT, nuclei)

# nuclei per subclass, one column per subclass -----------------------------
# subclasses are kept in LoadSubclassesToUse() order; replicates with zero
# nuclei in a subclass get 0 rather than NA.
df_subclass_counts <- df_meta |>
  as_tibble() |>
  filter(subclass_name %in% subclass_list) |>
  count(sample, subclass_name) |>
  mutate(subclass_name = factor(subclass_name, levels = subclass_list)) |>
  pivot_wider(
    names_from  = subclass_name,
    values_from = n,
    values_fill = 0,
    names_expand = TRUE
  )

df_supp_item1 <- df_sample |>
  left_join(df_subclass_counts, by = 'sample') |>
  mutate(across(all_of(subclass_list), \(x) replace_na(x, 0)))

stopifnot(
  ncol(df_supp_item1) == 8 + length(subclass_list),
  nrow(df_supp_item1) == n_distinct(df_meta$sample)
)

print(df_supp_item1, n = 50)

dir.create('05-results/SuppItem1', recursive = TRUE, showWarnings = FALSE)
write_csv(df_supp_item1, '05-results/SuppItem1/sample_metadata_subclass_counts.csv')