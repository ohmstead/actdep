## ----Fig2F-G.R
# This script will plot the expression of marker genes for dentate supertypes 
# in different activity conditions.
source('03-scripts/R/seq_functions.R')

if (!exists('nuclei')) {nuclei <- LoadDataset('Dec2024')}

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
# print(glue("{all_marker_genes[which(!all_marker_genes %in% rownames(mat))]} not in object"))
all_marker_genes <- all_marker_genes[which(all_marker_genes %in% rownames(mat))]

mat <- mat[all_marker_genes,] |> t()
df <- data.frame(mat) |> 
  rownames_to_column(var = 'barcode') |> 
  as_tibble() |> 
  left_join(nuclei_DG@meta.data) |> 
  pivot_longer(cols = any_of(all_marker_genes), names_to = 'gene', values_to = 'expression') |>
  relocate(gene, expression, .after = barcode)

# critical cluster-assignment genes
critical_genes <- c('Cenpa', 'Egr2', 'Egr4', # for cluster 0508
                    'Bhlhe41', 'Lct', 'Npy') # for cluster 0509

                    # activity_condition
p1 <- df |> 
  filter(gene %in% critical_genes) |>
  mutate(gene = factor(gene, levels = critical_genes)) |> 
  filter(str_detect(activity_condition, '^KA')) |>
ggplot() +
  aes(x = activity_condition, y = expression, fill = activity_condition) +
  geom_jitter(size = 3, width = 0.2, height = 0.1, shape = 21, alpha = 1, set.seed(17)) +
  # geom_jitter(width = 0.4, height = 0.1, shape = 21, size = 5, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, fill = 'gray80') +
  facet_wrap(~gene, scales = 'free_y') +
  scale_fill_manual(values = activity_colors) +
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 20, face = 'italic'),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
  )
print(p1)

if (exists('SAVE_PLOTS') & SAVE_PLOTS==TRUE) {
  ggsave(filename = 'critical__genes_activity.png', 
        path = '05-results/Figure2/raw_R_plots', 
        width = 8, height = 4, dpi = 900)
  ggsave(filename = 'critical_genes_activity.svg', 
        plot = p1 + LoadBarebonesTheme(ticks = 'y'),
        path = '05-results/Figure2/raw_R_plots', 
        width = 8, height = 4, dpi = 900)
}

# supertype
p2 <- df |> 
  filter(gene %in% critical_genes) |>
  mutate(gene = factor(gene, levels = critical_genes)) |> 
  filter(str_detect(cluster_name, '0508|0509')) |>
ggplot() +
  aes(x = cluster_name, y = expression, fill = cluster_name) +
  geom_jitter(size = 3, width = 0.4, height = 0.1, shape = 21, alpha = 1, set.seed(17)) +
  # geom_jitter(width = 0.4, height = 0.1, shape = 21, size = 5, set.seed(17)) +
  geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = NA, fill = 'gray80') +
  facet_wrap(~gene, scales = 'free_y') +
  scale_fill_manual(values = cluster_colors) +
  theme(
    legend.position = 'none',
    strip.text = element_text(size = 20, face = 'italic'),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
  )
print(p2)

if (exists('SAVE_PLOTS') & SAVE_PLOTS==TRUE) {
ggsave(filename = 'critical__genes_cluster.png', 
       path = '05-results/Figure2/raw_R_plots', 
       width = 8, height = 4, dpi = 900)
ggsave(filename = 'critical_genes_cluster.svg', 
       plot = p2 + LoadBarebonesTheme(ticks = 'y'),
       path = '05-results/Figure2/raw_R_plots', 
       width = 8, height = 4)
}
##----