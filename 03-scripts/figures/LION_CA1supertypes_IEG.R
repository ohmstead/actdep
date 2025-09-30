# This script plots percent activated cells in each CA1 supertype by activity condition.
# It also plots the distribution of active cells along anatomical axes.
source("03-scripts/R/seq_functions.R")

library(ggridges)
library(patchwork)
library(ggh4x)
library(gt)

condition_colors <- LoadActivityColors("May2024")
subclass_colors <- LoadAllenColors("subclass")
supertype_colors <- LoadAllenColors("supertype")
activity_colors <- LoadActivityColors()

nuclei <- LoadDataset("Dec2024")

# merge MERFISH meta ----------------------------------------
allen_taxonomy <- read_excel("02-data/published_data/Yao2023/allen_taxonomy_metadata.xlsx")
allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv")
meta_merfish <- read_csv("02-data/published_data/Zhang2023/cell_metadata.csv")

# IMPORTANT:
# "cluster_alias" col in merfish meta is equal to "cl" col in the allen_taxonomy
# "cl" and "cluster_id" are NOT the same thing in the allen_taxonomy. Use "cl" for joins!
allen_taxonomy <- allen_taxonomy |> 
  mutate(cl = as.numeric(cl))
ca1_merfish <- meta_merfish |> 
  left_join(allen_taxonomy, by = c("cluster_alias" = "cl")) |> 
  filter(low_quality_mapping == FALSE) |> 
  filter(subclass_id_label == '016 CA1-ProS Glut') |>
  filter(x < 9 | x > 2.5) |>      # anything outside this range is mis-classified
  filter(y < 8.1 | y > 2.9 ) |>   # anything outside is mis-classified
  filter(z < 8 & z > 4)           # anything outside is mis-classified

gene_list <- LoadGeneList("IEG")


# find active cells ------------------------------------------------
outputs <- FindActiveCells(nuclei, gene_list = gene_list, gene_threshold = 3)
df_active_cells_IEGs <- outputs$df_active_cells |>
  pivot_longer(cols = c(num_upregd_genes:last_col(), -num_upregd_genes), names_to = "gene", values_to = "expression") |> 
  left_join(nuclei@meta.data, by = c("cell" = "barcode", "activity_condition" = "activity_condition")) |> 
  filter(subclass_name == '016 CA1-ProS Glut')

df_supertype_active  <- df_active_cells_IEGs |>
  group_by(supertype_name, activity_condition) |> 
  summarise(fraction_active = sum(active_binary) / n())

df_activity_plotting_IEGs <- df_active_cells_IEGs |> 
  left_join(df_supertype_active, by = c("supertype_name", "activity_condition")) |> 
  summarize(
    .by = cell,
    supertype_name = dplyr::first(supertype_name), 
    num_upregd_genes = dplyr::first(num_upregd_genes), 
    activity_condition = dplyr::first(activity_condition)
  )


# p_supertypes_IEG  ----------------------------------------
strip <- strip_themed(background_x = elem_list_rect(fill = supertype_colors[c(101, 107, 102, 108, 106, 97)]))
strip <- strip_themed(background_x = elem_list_rect(fill = activity_colors))

p_supertypes_IEG <- ggplot(df_activity_plotting_IEGs) + 
  aes(y = num_upregd_genes, x = supertype_name, color = supertype_name) +
  geom_jitter( size = 0.05, alpha = 0.3, width = 0.3, height = 0.25, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, fill = 'gray80', color = 'black') +
  geom_hline(yintercept = 2.5, linetype = 'dotted') +
  scale_color_manual(values = supertype_colors) +
  scale_y_continuous(breaks = c(0, 3, 5, 10, 15)) +
  theme_void() +
  theme(
    legend.position = 'none',
    strip.text = element_blank(),
  ) +
  facet_wrap2(~activity_condition, strip = strip, nrow = 1)

print(p_supertypes_IEG)


# table ------------------------------------------
# remap supertype names for the table
supertype_map <- c(
  "0069 CA1-ProS Glut_1" = "CA1 69-1",
  "0070 CA1-ProS Glut_2" = "CA1 70-2",
  "0071 CA1-ProS Glut_3" = "CA1 71-3",
  "0072 CA1-ProS Glut_4" = "CA1 72-4",
  "0073 CA1-ProS Glut_5" = "CA1 73-5",
  "0074 CA1-ProS Glut_6" = "CA1 74-6"
)
ca1_supertype_colors <- supertype_colors[sort(unique(df_activity_plotting_IEGs$supertype_name))]
names(ca1_supertype_colors) <- recode(names(ca1_supertype_colors), !!!supertype_map)

gt_summary <- df_activity_plotting_IEGs |> 
  mutate(supertype_name = recode(supertype_name, !!!supertype_map)) |>
  add_count(supertype_name, name = "Total n") |> 
  group_by(supertype_name, activity_condition, `Total n`) |> 
  summarize(pct_active = sum(num_upregd_genes > 2.5) / n(), .groups = "drop") |> 
  pivot_wider(names_from = activity_condition, values_from = pct_active) |> 
  ungroup() |> 
gt() |>  # Do NOT set rowname_col here
  tab_header(title = "Percent of active CA1 cells by Supertype") |> 
  fmt_percent(columns = -`Total n`, decimals = 0) |> 
  fmt_number(columns = `Total n`, use_seps = TRUE, decimals = 0) |>
  opt_table_font(font = "Helvetica") |> 
  cols_label(supertype_name = '') |> 
  cols_align(align = "left", columns = everything()) |>
  # cols_width(2 ~ px(60)) |>
  # cols_width(3 ~ px(70)) |>
  # cols_width(4 ~ px(65)) |>
  # cols_width(5 ~ px(70)) |>
  # cols_width(6 ~ px(70)) |>
  data_color(
    columns = "supertype_name",
    fn = scales::col_factor(
      palette = ca1_supertype_colors,
      domain  = names(ca1_supertype_colors)
    ),
    apply_to = "text"
  ) |> 
  tab_style(
    style    = cell_text(size = px(16)),
    locations = cells_body()
  )

print(gt_summary)
gtsave(gt_summary, "05-results/LION/raw_R_plots/CA1_supertypes_activation_table.png")
gtsave(gt_summary, "05-results/LION/raw_R_plots/CA1_supertypes_activation_table.pdf")


# p_CA1_supertypes_APaxis_IEG ------------------------------------------------
p_CA1_supertypes_APaxis_IEG <- ca1_merfish |> 
  group_by(z, supertype_id_label) |>  
  summarise(n = n()) |> 
  mutate(composition = n / sum(n)) |> 
  arrange(z, supertype_id_label) |> 
  left_join(df_supertype_active, by = c("supertype_id_label" = "supertype_name")) |> 
  mutate(height = composition * fraction_active) |> 
ggplot() +
  geom_area(aes(x = z, y = height, fill = supertype_id_label), position = 'stack') +
  geom_vline(xintercept = seq(7.5, 4.2, -0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = supertype_colors) +
  scale_x_reverse() +
  theme_void() +
  theme(
    # strip.text = element_text(size = 20, color = 'white'),
    strip.text = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none',
    legend.text = element_text(size = 20)
  ) +
  facet_wrap2(~activity_condition, strip = strip, nrow = 1)
print(p_CA1_supertypes_APaxis_IEG)


# p_CA1_supertypes_DVaxis_IEG ------------------------------------------------
p_CA1_supertypes_DVaxis_IEG <- ca1_merfish |> 
  mutate(ycut = cut(y, breaks = seq(2.9, 8.1, 0.2))) |> 
  mutate(yy = as.numeric(substr(as.character(ycut), 2, 4))) |>   # change to number
  group_by(yy, supertype_id_label) |> 
  summarise(n = n()) |>
  mutate(composition = n / sum(n)) |> 
  arrange(yy, supertype_id_label) |> 
  left_join(df_supertype_active, by = c("supertype_id_label" = "supertype_name")) |> 
  mutate(height = composition * fraction_active) |>
ggplot() +
  geom_area(aes(x = yy, y = height, fill = supertype_id_label), position = 'stack') +
  geom_vline(xintercept = seq(3, 8, 0.2), color = 'white', linetype = 2, alpha = 0.3) +
  scale_fill_manual(values = supertype_colors) +
  coord_flip() +
  scale_x_reverse() +
  theme_void() +
  theme(
    # strip.text = element_text(size = 20, color = 'white'),
    strip.text = element_blank(),
    plot.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none',
    legend.text = element_text(size = 20)
  ) +
  facet_wrap2(~activity_condition, strip = strip, nrow = 1)
print(p_CA1_supertypes_DVaxis_IEG)


if (SAVE_PLOTS == TRUE) {
  save_plots <- "05-results/LION/raw_R_plots"
  # png
  ggsave(plot = p_supertypes_IEG, 
         path = save_plots, 
         filename = "CA1_supertypes_IEG.png",
         device = png, width = 16, height = 5, dpi = 300)
  ggsave(plot = p_IEG_spatial_AP,
         path = save_plots,
         filename = "CA1_supertypes_APaxis_IEG.png",
         device = png, width = 9, height = 1.5, dpi = 300)
  ggsave(plot = p_IEG_spatial_DV,
         path = save_plots,
         filename = "CA1_supertypes_DVaxis_IEG.png",
         device = png, width = 16, height = 4, dpi = 300)
  # svg
  ggsave(plot = p_supertypes_IEG, 
         path = save_plots, 
         filename = "CA1_supertypes_IEG.svg",
         width = 9, height = 3, units = 'in')
  ggsave(plot = p_CA1_supertypes_APaxis_IEG,
         path = save_plots,
         filename = "CA1_supertypes_APaxis_IEG.svg",
         width = 9, height = 1.5, units = 'in')
  ggsave(plot = p_CA1_supertypes_DVaxis_IEG,
         path = save_plots,
         filename = "CA1_supertypes_DVaxis_IEG.svg",
         width = 9, height = 1, units = 'in')
}
