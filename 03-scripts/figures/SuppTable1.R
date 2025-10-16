# make a gt table of sample metadata
library(gt)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset('Dec2024')
activity_colors <- LoadActivityColors()
zt_colors <- LoadZTColors()
sex_colors <- LoadSexColors()

# load and clean sample meta; count nuclei
df_metadata_csv <- read_csv("01-documentation/sample_experimental_metadata_Dec2024.csv")
df_sample <- nuclei@meta.data |> 
  as_tibble() |> 
  select(sample, age, sex, hemisphere_seq, ZT, activity_condition) |> 
  left_join(df_metadata_csv, by = 'sample', suffix = c('.seurat', '')) |> 
  group_by(sample) |> 
  summarize(
    sample = first(sample),
    age = first(age),
    sex = first(sex),
    hemisphere = first(hemisphere_seq),
    ZT = first(ZT),
    activity_condition = first(activity_condition),
    nuclei = n()
  ) |> 
  arrange(activity_condition, sample) |> 
  print(n=50)


# condition x ZT ----------------------------------------
table_conditions <- df_sample |> 
  # mutate(sample_num = 1:length(sample)) |>
  group_by(activity_condition, ZT) |> 
  summarize(
    age_sd = sd(age),
    age = mean(age),
    Replicates = n(),
    male = sum(sex=='M'),
    female = sum(sex=='F'),
    R = sum(hemisphere=='R'),
    L = sum(hemisphere=='L'),
    nuclei_total = sum(nuclei),
    nuclei_mean = mean(nuclei),
    nuclei_sd = sd(nuclei),
    .groups = 'drop'
    ) |> 
  arrange(activity_condition) |> 
  select(-nuclei_sd) |>
gt(groupname_col = 'activity_condition', rowname_col = 'ZT', row_group_as_column = T, rownames_to_stub = F) |> 
  tab_header(title = 'Factor design in ACT-DEPP dataset') |> 
  fmt_number(columns = c(nuclei_total, nuclei_mean), use_seps = T, decimals = 0) |> 
  fmt_number(columns = c(age, age_sd), use_seps = T, decimals = 1) |>
  tab_spanner(label = 'Sex', columns = c(male, female)) |>
  tab_spanner(label = 'Hemisphere', columns = c(R, L)) |>
  tab_spanner(label = 'Nuclei', columns = c(nuclei_total, nuclei_mean)) |> 
  cols_merge_uncert(col_val = age, col_uncert = age_sd) |> 
  cols_align(align = 'left', columns = ZT) |> 
  sub_missing(columns = ZT, missing_text = "---") |> 
  cols_label(
    male = 'M',
    female = 'F',
    age = 'Age ± sd',
    nuclei_mean = 'Average',
    nuclei_total = 'Total',
  ) |>
  cols_width(
    ZT ~ px(60),
    c(male, female) ~ px(50),
    nuclei_total ~ px(80),
    nuclei_mean ~ px(80),
    ) |> 
  summary_rows(
    groups = c('SE','EE30m', 'EE6h'),
    columns = c(Replicates, male, female, R, L, nuclei_total, nuclei_mean),
    fns = list(label = md(glue('_n~total~_'))) ~ sum(.),
    fmt = ~fmt_number(., use_seps = T, decimals = 0),
    missing_text = ''
  ) |> 
  grand_summary_rows(
    columns = c(Replicates, male, female, R, L, nuclei_total),
    fns = list(label = md(glue('_n<sub>ACT-DEPP</sub>_'))) ~ sum(.),
    fmt = ~fmt_number(., use_seps = T, decimals = 0),
    missing_text = ''
  ) |> 
  cols_align(
    align = "center",
    columns = c(age, Replicates, male, female, R, L, nuclei_total, nuclei_mean)
  ) |> 
  cols_align(
    align = "right",
    columns = c(nuclei_total, nuclei_mean)
  ) |> 
  tab_style(
    style = cell_text(align = "center"),
    locations = cells_column_spanners(
      spanners = c("Sex", "Hemisphere", "Nuclei")
    )
  )|> 
  opt_table_font(font = 'Helvetica')
  
print(table_conditions)
gtsave(path = '05-results/SuppTable1/raw_R_plots/', 
       filename = 'condition_table.png',
       table_conditions)

  

# all_samples ----------------------------------------
tbl_sample <- df_sample |> 
  mutate(sample_num = 1:length(sample)) |>
  select(-c(sample_num, sample)) |> 
  mutate(nuclei = as.numeric(nuclei)) |> 
gt(groupname_col = 'activity_condition', row_group_as_column = T, rownames_to_stub = T) |> 
  tab_header(title = 'Biological replicates in dataset') |> 
  fmt_number(columns = c(nuclei), use_seps = T, decimals = 0) |>
  sub_missing(columns = ZT, missing_text = "---") |> 
  summary_rows(
    groups = everything(),
    columns = nuclei,
    fns = list(label = md('_n~replicates~_')) ~n(),
    missing_text = '',
  ) |> 
  summary_rows(
    groups = everything(),
    columns = nuclei,
    fns = list(label = md('_n~nuclei~_')) ~sum(.),
    fmt = ~fmt_number(., use_seps = T, decimals = 0),
    missing_text = ''
  ) |> 
  grand_summary_rows(
    columns = nuclei,
    fns = list(label = md('_N<sub>total replicates</sub>_')) ~n(),
    missing_text = ''
  ) |> 
  grand_summary_rows(
    columns = nuclei,
    fns = list(label = md('_N<sub>total nuclei</sub>_')) ~sum(.),
    fmt = ~fmt_number(., use_seps = T, decimals = 0),
    missing_text = ''
  ) |> 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = activity_condition)
  )
print(tbl_sample)

gtsave(path = '05-results/SuppTable1/raw_R_plots/', 
       filename = 'all_samples.png',
       tbl_sample)
