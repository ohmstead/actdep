# load libs & data -----
source("03-scripts/R/seq_functions.R")

library(dplyr)
library(readr)
library(readxl)
library(ggridges)
library(patchwork)
library(ggplot2)
library(ggh4x)
library(gt)

condition_colors <- LoadActivityColors("May2024")
subclass_colors <- LoadAllenColors("subclass")
supertype_colors <- LoadAllenColors("supertype")
activity_colors <- LoadActivityColors()

nuclei <- LoadDataset("Dec2024")

# load MERFISH and metadata ----
allen_taxonomy <- read_excel("02-data/published_data/Yao2023/allen_taxonomy_metadata.xlsx")
allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv")
meta_merfish <- read_csv("02-data/published_data/Zhang2023/cell_metadata.csv")


# merge MERFISH and metadata ----
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


# find IEG activated cells ----
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


# plot IEG activation ---- 
strip <- strip_themed(background_x = elem_list_rect(fill = supertype_colors[c(101, 107, 102, 108, 106, 97)]))
strip <- strip_themed(background_x = elem_list_rect(fill = activity_colors))

p_supertypes_IEG <- ggplot(df_activity_plotting_IEGs) + 
  aes(y = num_upregd_genes, x = supertype_name, fill = supertype_name) +
  geom_jitter(width = 0.3, height = 0.25, shape = 21, alpha = 0.3, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, fill = 'gray80') +
  geom_hline(yintercept = 2.5, linetype = 'dotted') +
  scale_fill_manual(values = supertype_colors) +
  scale_y_continuous(breaks = c(0, 3, 5, 10, 15)) +
  theme(
    plot.title = element_text(size = 30),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none',
    # strip.text = element_text(size = 15, color = 'white'),
    strip.text = element_blank(),
    legend.text = element_text(size = 15)
  ) +
  facet_wrap2(~activity_condition, strip = strip, nrow = 1)

print(p_supertypes_IEG)

# print table
ca1_supertype_colors <- supertype_colors[sort(unique(df_activity_plotting_IEGs$supertype_name))]
gt_summary <- df_activity_plotting_IEGs |> 
  group_by(supertype_name, activity_condition) |> 
  summarize(pct_active = sum(num_upregd_genes > 2.5) / n()) |>
  pivot_wider(names_from = activity_condition, values_from = pct_active) |> 
  ungroup() |> 
gt(rowname_col = 'supertype_name') |> 
  tab_header(title = "Percent of cells active",
             subtitle = 'CA1 supertypes') |> 
  fmt_percent(decimals = 1) |> 
  opt_table_font(font = "Aptos") |> 
  tab_style(
    style = cell_fill(color = subclass_colors[match(.$supertype_name, names(subclass_colors))]), 
    locations = cells_stub(rows = everything())  # Apply color to rowname column
  )
print(gt_summary)
gtsave(gt_summary, "05-results/figure2/raw_R_plots/CA1_supertypes_IEG_summary.pdf")

# plot activation along anatomical axes ----
# A/P
p_IEG_spatial_AP <- ca1_merfish |> 
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
print(p_IEG_spatial_AP)

# D/V
p_IEG_spatial_DV <- ca1_merfish |> 
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
print(p_IEG_spatial_DV)


if (SAVE_PLOTS == TRUE) {
  ggsave(plot = p_supertypes_IEG, 
         path = "05-results/figure2/raw_R_plots", 
         filename = "CA1_supertypes_IEG.png",
         device = png, width = 16, height = 5, dpi = 300)
  ggsave(plot = p_IEG_spatial_AP,
         path = "05-results/figure2/raw_R_plots",
         filename = "CA1_supertypes_APaxis_IEG.png",
         device = png, width = 16, height = 5, dpi = 300)
  ggsave(plot = p_IEG_spatial_DV,
         path = "05-results/figure2/raw_R_plots",
         filename = "CA1_supertypes_DVaxis_IEG.png",
         device = png, width = 16, height = 4, dpi = 300)
}
