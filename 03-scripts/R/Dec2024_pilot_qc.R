library(tidyverse)
library(plotly)

# for each sublibrary ----
# for each sample, plot the stability of number of cells per sublibrary. requires a df with info about samples per sublibrary
parent_dir <- "02-data/raw_data/sublibrary_reports"

# go into each subdir, get the name
sublibrary_dirs <- list.dirs(parent_dir, full.names = FALSE, recursive = FALSE) |> 
  map_chr(~str_remove(.x, parent_dir)) |> 
  str_remove("^/") |> 
  str_remove("/$")

# for each sublibrary_dir, read in "agg_samp_ana_summary.csv" and add a column for sublibrary
per_sublibrary <- map_dfr(sublibrary_dirs, ~{
  sublibrary <- .x
  read_csv(file.path(parent_dir, sublibrary, "agg_samp_ana_summary.csv")) |> 
    pivot_longer(cols = -statistic, names_to = 'sample_name') |> 
    mutate(condition = str_extract(sample_name, "^[^_]+")) |>
    mutate(sample = str_remove(sample_name, "^[^_]+_")) |> 
    relocate(condition, sample, .after = sample_name) |> 
    pivot_wider(names_from = statistic, values_from = value) |> 
    select(-starts_with("GRCm39")) |> 
    filter(condition != 'all-sample') |> 
    group_by(sample_name) |> 
    mutate(cells_per_1M_reads = number_of_cells / number_of_reads * 1e6) |>
    mutate(tscp_per_1M_reads = number_of_tscp / number_of_reads * 1e6) |> 
    mutate(sublibrary = str_extract(sublibrary, "JO\\d+")) |> 
    relocate(sublibrary) |> 
    ungroup()
}) |> print()

p <- ggplot(per_sublibrary) +
  aes(x = sublibrary, y = number_of_cells, fill = condition, z = sample_name) +
  geom_boxplot(aes(fill = 'white'), outlier.shape = NA) +
  # geom_jitter(shape = 21, size = 3, width = 0.2, alpha = 0.8) +
  geom_point(shape = 21, size = 3, alpha = 0.8) +
  geom_line(aes(group = sample_name), alpha = 0.1) +
  geom_smooth(aes(group = condition, color = condition), method = 'loess', se = FALSE, linetype = 'dashed') +
  labs(title = 'Number of cells') +
  theme(title = element_text(size = 30), legend.position = 'none')
ggplotly(p)
# 
# p <- ggplot(per_sublibrary) +
#   aes(x = sublibrary, y = number_of_reads, fill = sublibrary, z = sample_name) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(shape = 21, size = 3, width = 0.2, alpha = 0.3) +
#   labs(title = 'Number of cells') +
#   theme(title = element_text(size = 30), legend.position = 'none')
# ggplotly(p)
# 
# p <- ggplot(per_sublibrary) +
#   aes(x = sublibrary, y = cells_per_1M_reads, fill = condition, z = sample_name) +
#   geom_boxplot(outlier.shape = NA) +
#   geom_jitter(shape = 21, size = 3, width = 0.2, alpha = 0.3) +
#   labs(title = 'Number of cells') +
#   theme(title = element_text(size = 30), legend.position = 'none')
# ggplotly(p)
