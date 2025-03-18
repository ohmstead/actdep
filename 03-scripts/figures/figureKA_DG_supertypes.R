# This script will plot the expression of marker genes for dentate supertypes 
# in different activity conditions.

library(ggplot2)
library(ggvenn)
library(patchwork)
library(readr)
library(dplyr)
library(readxl)
library(glue)

library(Seurat)

source('03-scripts/R/seq_functions.R')

nuclei <- LoadDataset('Dec2024')
activity_colors  <- LoadActivityColors()
subclass_colors  <- LoadAllenColors('subclass')
supertype_colors <- LoadAllenColors('supertype')
cluster_colors   <- LoadAllenColors('cluster')


# get marker genes ------------------------------------------------
supertype_sheet <- read_excel("02-data/published_data/Yao2023/allen_taxonomy_metadata.xlsx",
                  sheet = "supertype_annotation")
cluster_sheet <- read_excel("02-data/published_data/Yao2023/allen_taxonomy_metadata.xlsx",
                  sheet = "cluster_annotation")
supertype_genes_DG <- supertype_sheet |> 
  slice(136:139) |> 
  select(contains('markers.combo')) |>
  pivot_longer(cols = everything(), values_to = 'gene', names_to = NULL) |> 
  separate_rows(gene, sep = ",") |> 
  distinct(gene) |> 
  pull(gene)

cluster_genes_DG <- cluster_sheet |> 
  slice(502:510) |> 
  select(contains('markers.combo')) |>
  pivot_longer(cols = everything(), values_to = 'gene', names_to = NULL) |> 
  separate_rows(gene, sep = ",") |> 
  distinct(gene) |> 
  pull(gene)



# df expression ---------------------------------------------
nuclei_DG <- nuclei |> subset(subclass_name == '037 DG Glut')

mat <- GetAssayData(nuclei_DG, assay = 'SCT', layer = 'data')
all_marker_genes <- unique(c(supertype_genes_DG, cluster_genes_DG))
print(glue("{all_marker_genes[which(!all_marker_genes %in% rownames(mat))]} not in object"))
all_marker_genes <- all_marker_genes[which(all_marker_genes %in% rownames(mat))]

mat <- mat[all_marker_genes,] |> t()
df <- data.frame(mat) |> 
  rownames_to_column(var = 'barcode') |> 
  as_tibble() |> 
  left_join(nuclei_DG@meta.data) |> 
  pivot_longer(cols = any_of(all_marker_genes), names_to = 'gene', values_to = 'expression') |>
  relocate(gene, expression, .after = barcode) |> 
  print()


# plot -------------------------------------------------------
# supertype genes
df |> 
  filter(gene %in% supertype_genes_DG) |>
ggplot() +
  aes(x = supertype_name, y = expression, fill = activity_condition) +
  geom_violin() +
  facet_wrap(~gene, scales = 'free_y') +
  scale_fill_manual(values = activity_colors)

# cluster genes: activity
df |> 
  filter(gene %in% cluster_genes_DG) |>
ggplot() +
  aes(x = activity_condition, y = expression, fill = activity_condition) +
  geom_boxplot(outliers = F) +
  # geom_violin() +
  facet_wrap(~gene, scales = 'free_y') +
  scale_fill_manual(values = activity_colors)

# cluster genes: supertype
df |> 
  filter(gene %in% cluster_genes_DG) |>
ggplot() +
  aes(x = supertype_name, y = expression, fill = supertype_name) +
  geom_boxplot(outliers = F) +
  facet_wrap(~gene, scales = 'free_y') +
  scale_fill_manual(values = supertype_colors)

DoHeatmap(nuclei_DG, 
          features = cluster_genes_DG, 
          group.by = 'cluster_name', 
          slot = 'data',
          group.colors = cluster_colors[sort(unique(df$cluster_name))],
          group.bar.height= 0.05) +
  scale_fill_viridis_c(limits = c(0,2), oob = scales::squish)


# critical cluster-assignment genes
critical_genes <- c('Cenpa', 'Egr2', 'Egr4', # for cluster 0508
                    'Bhlhe41', 'Lct', 'Npy') # for cluster 0509
# activity_condition
df |> 
  filter(gene %in% critical_genes) |>
  mutate(gene = factor(gene, levels = critical_genes)) |> 
ggplot() +
  aes(x = activity_condition, y = expression, fill = activity_condition) +
  # geom_jitter(width = 0.4, height = 0.1, shape = 21, alpha = 0.3, set.seed(17)) +
  geom_jitter(width = 0.4, height = 0.1, shape = 21, size = 5, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, fill = 'gray80') +
  facet_wrap(~gene, scales = 'free_y') +
  scale_fill_manual(values = activity_colors) +
  theme(
    # legend.position = 'none',
    strip.text = element_text(size = 20),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
  )
si(1000, 500, format = 'png')

# supertype
df |> 
  filter(gene %in% critical_genes) |>
  mutate(gene = factor(gene, levels = critical_genes)) |> 
  ggplot() +
  aes(x = cluster_name, y = expression, fill = cluster_name) +
  geom_jitter(width = 0.4, height = 0.1, shape = 21, alpha = 0.3, set.seed(17)) +
  # geom_jitter(width = 0.4, height = 0.1, shape = 21, size = 5, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, fill = 'gray80') +
  facet_wrap(~gene, scales = 'free_y') +
  scale_fill_manual(values = cluster_colors) +
  theme(
    # legend.position = 'none',
    strip.text = element_text(size = 20),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
  )
si(1000, 500, format = 'png')




# 0503 DG Glut_1  0504 DG Glut_1  0505 DG Glut_2  0506 DG Glut_2  0507 DG Glut_2  0508 DG Glut_3  0509 DG Glut_3  0510 DG Glut_4 
c("#5CAECC",      "#FF0700",      "#CC603D",      "#FFD826",      "#99176A",      "#2E9964",      "#CC3DA3",      "#8EFF4D" )
