# load libs & data -----
source("03-scripts/R/seq_functions.R")

library(dplyr)
library(readr)
library(readxl)
library(ggridges)
library(patchwork)
library(ggplot2)
library(ggh4x)
library(plotly)

condition_colors <- LoadActivityColors("May2024")
subclass_colors <- LoadAllenColors("subclass")
supertype_colors <- LoadAllenColors("supertype")

# nuclei <- LoadDataset("May2024", "combined")

# load MERFISH and metadata ----
allen_taxonomy <- read_excel("02-data/published_data/allen_taxonomy_metadata.xlsx")
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


# find IEG activated cells ----
gene_list <- LoadGeneList("IEG")
fxn_outputs <- FindActiveCells(nuclei, gene_list = gene_list, gene_threshold = 3)
df_active_cells_IEGs <- fxn_outputs$df_active_cells |>
  pivot_longer(cols = c(num_upregd_genes:last_col(), -num_upregd_genes), names_to = "gene", values_to = "expression") |> 
  mutate(is_expressing = expression > 0) |> 
  right_join(nuclei@meta.data, by = c("cell" = "barcode", "activity_condition" = "activity_condition")) |> 
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
p_IEG_active_supertype <- ggplot(df_activity_plotting_IEGs) + 
  aes(y = num_upregd_genes, x = supertype_name, fill = supertype_name) +
  geom_jitter(width = 0.3, height = 0.25, shape = 21, alpha = 0.3, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 2.5, linetype = 'dotted') +
  scale_fill_manual(values = supertype_colors) +
  scale_y_continuous(breaks = c(0, 3, 5, 10, 15)) +
  theme(
    plot.title = element_text(size = 30),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none',
    strip.text = element_text(size = 15, color = 'white'),
    legend.text = element_text(size = 15)
  ) +
  facet_wrap2(~activity_condition, strip = strip, nrow = 1)
  
print(p_IEG_active_supertype)

nuclei@meta.data |> 
  filter(subclass_name %in% subclass_list$subclass_name) |> 
  summarize(n = n(), .by = c('subclass_name', 'activity_condition')) |> 
  arrange(subclass_name, n)


# plot activation along anatomical axes ----
p_IEG_spatial <- ca1_merfish |> 
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
    strip.text = element_text(size = 20, color = 'white'),
    plot.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = 'none',
    legend.text = element_text(size = 20)
  ) +
  facet_wrap2(~activity_condition, strip = strip, nrow = 1)

print(p_IEG_spatial)


if (SAVE_PLOTS == TRUE) {
  ggsave(plot = p_IEG_active_supertype, 
         path = "05-results/figure2/raw_R_plots", 
         filename = "IEG_activation_by_supertype.png",
         device = png, width = 16, height = 5, dpi = 900)
  ggsave(plot = p_IEG_spatial,
         path = "05-results/figure2/raw_R_plots",
         filename = "IEG_activation_by_APaxis.png",
         device = png, width = 16, height = 5, dpi = 900)
  
}


# plot co-expression by supertype ----
EE30m_6iegs <- df_active_cells |> 
  filter(activity_condition == 'EE30m') |>
  # filter(num_upregd_genes == 6) |> 
  filter(supertype_name == '0072 CA1-ProS Glut_4')

df_active_cells |> 
  group_by(cell) |> 
  summarize(
    type = first(supertype_name),
    n.genes = first(num_upregd_genes),
    cond = first(activity_condition)
  ) |> 
ggplot() +
  aes(x = n.genes, fill = cond) +
  geom_histogram() +
  facet_grid(vars(type), vars(cond), scales = 'free_y')
