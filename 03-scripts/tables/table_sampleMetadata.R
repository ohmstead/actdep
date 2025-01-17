# This script creates a sample metadata table for the supplement
library(gt)
library(tidyverse)

source("03-scripts/R/seq_functions.R")

# load in data
activity_colors <- LoadActivityColors()
zt_colors <- LoadZTColors()
sex_colors <- LoadSexColors()

meta <- read_csv("01-documentation/sample_experimental_metadata_Dec2024.csv")

# mutate new columns ----
meta <- meta |> 
  mutate(ZT = str_extract(sample, "^[^_]+")) |>
  mutate(seizure = ifelse(str_detect(ZT, "KA"), TRUE, FALSE)) |>
  mutate(ZT = ifelse(seizure == T, NA, ZT)) |> 
  mutate(age = paste0("P", age))

# make col: activity_condition (e.g. SE, EE6h, KA30m)
meta <- meta |>
  mutate(activity_condition = str_sub(sample, 1, -2)) |>  # remove replicate numbers
  mutate(activity_condition = ifelse(
    str_detect(str_sub(activity_condition, 1, 2), "KA"),         # if KA
    str_remove(activity_condition, '_'),                         # then make it KA30m or KA6h
    str_remove(str_extract(activity_condition, "(?<=_).*"), '_') # else only include text after first underscore
  )) |> 
  mutate(activity_condition = ifelse(         # add "EE_" prefix to 30m and 6h
    str_detect(activity_condition, "SE"), 
    "SE", 
    ifelse(
      !str_detect(activity_condition, "KA"), 
      paste0("EE", activity_condition), 
      activity_condition
    )))

# create levels to factors
meta$ZT <- factor(meta$ZT, levels = c('ZT0', 'ZT4', 'ZT12', 'ZT16'))
meta$activity_condition <- factor(meta$activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))

# arrange rows for table order
meta <- meta |> 
  relocate(activity_condition, ZT, sex, age, hemisphere_seq, sample) |> 
  arrange(activity_condition, ZT) |> 
  rename(
    "Sex" = sex,
    "Age" = age,
    "Hemisphere" = hemisphere_seq
  )


# make gt table ----
meta |> 
  select(-sample) |> 
  select(-seizure) |> 
# gt(rowname_col = 'activity_condition') |> 
gt() |> 
  tab_header(title = "Sample metadata") |> 
  cols_label(activity_condition = "") |>
  tab_row_group(
    label = 'Physiological',
    rows = 1:42
  ) |> 
  tab_row_group(
    label = 'Seizure',
    rows = 43:48
  ) |> 
  row_group_order(
    groups = c('Physiological', 'Seizure')
  ) |> 
  data_color(
    columns = Sex,
    palette = sex_colors,
    domain = c('M', 'F'),
    ordered = T
    # fn = scales::col_factor(sex_colors, domain = names(sex_colors))
  ) |> 
  data_color(
    columns = activity_condition,
    palette = activity_colors,
    domain = c("SE", "EE30m", "EE6h", "KA30m", "KA1h", "KA6h"),
    ordered = T
  ) |> 
  data_color(
    columns = ZT,
    palette = zt_colors,
    domain = c('ZT0', 'ZT4', 'ZT12', 'ZT16'),
    ordered = T
  ) |> 
  gtsave("01-documentation/table_sampleMetadata.docx")