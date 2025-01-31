library(tidyverse)
library(ggsankey)
library(reticulate)
library(scCustomize)

reticulate::use_condaenv("sc")
source("03-scripts/R/seq_functions.R")
plot_save_dir <- '05-results/figure1/raw_R_plots'

# load in seurat
nuclei_with <- LoadDataset("Dec2024")

# 1. Run Dec2024 with ARGs through MapMyCells.
      # already done

# 2. Remove Tyssowski ARGs from Seurat object.
genes_to_remove <- LoadGeneList("tyssowski")
features_with <- Features(nuclei_with)
features_without <- features_with[!features_with %in% genes_to_remove]

# use DietSeurat to get just the RNA assay, which we will then subset
DefaultAssay(nuclei_with) <- 'RNA'
nuclei_without <- DietSeurat(nuclei_with, 
                             assays = 'RNA', 
                             layers = 'counts', 
                             features = features_without)

# Run Dec2024 without ARGs through MapMyCells ---- 
working_save_dir <- "04-analysis/cluster_reassignment"
as.anndata(x = nuclei_without, file_path = working_save_dir, file_name = "Dec2024_withoutARGs", main_layer = "counts", other_layers = NULL)


# merge with and without metadata ----
meta_mapmycells <- read_csv("04-analysis/cluster_reassignment/Dec2024_withoutARGs_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1738267008578/Dec2024_withoutARGs_10xWholeMouseBrain(CCN20230722)_HierarchicalMapping_UTC_1738267008578.csv", skip = 4)
meta_merge <- nuclei_with@meta.data |> 
      left_join(meta_mapmycells,
                by = c('barcode' = 'cell_id'),
                suffix = c('.with', '.without'),
                relationship = 'one-to-one',
                unmatched = 'error'
                )

# make some new vars
subclasses_to_use <- LoadSubclassesToUse(nuclei_with)
meta <- meta_merge |> 
      mutate(class_bootstrapping_probability_diff = class_bootstrapping_probability.with - class_bootstrapping_probability.without) |> 
      mutate(subclass_bootstrapping_probability_diff = subclass_bootstrapping_probability.with - subclass_bootstrapping_probability.without) |> 
      mutate(supertype_bootstrapping_probability_diff = supertype_bootstrapping_probability.with - supertype_bootstrapping_probability.without) |> 
      mutate(cluster_bootstrapping_probability_diff = cluster_bootstrapping_probability.with - cluster_bootstrapping_probability.without) |> 
      filter(subclass_name.with %in% subclasses_to_use$subclass_name) |> 
      select(-starts_with('emptyDrops_'))


# plot difference in mapping probability with/without ARGs
class_colors = LoadAllenColors('class')
subclass_colors = LoadAllenColors('subclass')
supertype_colors = LoadAllenColors('supertype')
cluster_colors = LoadAllenColors('cluster')

# class
meta |> 
ggplot() +
      aes(x = class_name.with, y = class_bootstrapping_probability_diff, fill = class_name.with) +
      geom_jitter(shape = 21) +
      geom_boxplot(alpha = 0.5, outliers = F) +
      scale_fill_manual(values = class_colors)

# subclass
meta |> 
ggplot() +
      aes(x = subclass_name.with, y = subclass_bootstrapping_probability_diff, fill = subclass_name.with) +
      geom_jitter(shape = 21) +
      geom_boxplot(alpha = 0.5, outliers = F) +
      scale_fill_manual(values = subclass_colors)

# supertype
meta |> 
      filter(subclass_name.with == '016 CA1-ProS Glut') |> 
ggplot() +
      aes(x = supertype_name.with, y = supertype_bootstrapping_probability_diff, fill = supertype_name.with) +
      geom_jitter(alpha = 0.3, shape = 21) +
      geom_boxplot(alpha = 0.5, outliers = F) +
      scale_fill_manual(values = supertype_colors)

# cluster
meta |> 
      filter(subclass_name.with == '016 CA1-ProS Glut') |> 
ggplot() +
      aes(x = cluster_name.with, y = cluster_bootstrapping_probability_diff, fill = cluster_name.with) +
      geom_jitter(shape = 21) +
      geom_boxplot(alpha = 0.5, outliers = F) +
      scale_fill_manual(values = cluster_colors)


# find difference in CA1 active cells ----
IEGs <- LoadGeneList()
outputs <- FindActiveCells(nuclei_with, gene_list = IEGs)
df_active_cells <- outputs$df_active_cells |> 
      select(cell, active_binary)
meta <- meta |> 
      left_join(df_active_cells,
                by = c('barcode' = 'cell'),
                suffix = c('.with', '.without'),
                unmatched = 'error'
                )
active_binary_colors <- setNames(c('red', 'black'), c(T, F))

# supertype
meta |> 
      filter(subclass_name.with == '016 CA1-ProS Glut') |> 
ggplot() +
      aes(x = supertype_name.with, y = supertype_bootstrapping_probability_diff, 
          fill = supertype_name.with, color = active_binary) +
      geom_jitter(alpha = 0.3) +
      geom_boxplot(alpha = 0.5, outliers = F) +
      scale_fill_manual(values = supertype_colors) +
      scale_color_manual(values = active_binary_colors)

# cluster
meta |> 
      filter(subclass_name.with == '016 CA1-ProS Glut') |> 
ggplot() +
      aes(x = cluster_name.with, y = cluster_bootstrapping_probability_diff, 
          fill = cluster_name.with, color = active_binary) +
      geom_jitter(alpha = 0.3) +
      geom_boxplot(alpha = 0.5, outliers = F) +
      scale_fill_manual(values = cluster_colors) +
      scale_color_manual(values = active_binary_colors)


# sankey diagrams of re-assignment ----
# class sankey
meta |> 
      make_long(class_name.with, class_name.without) |> 
ggplot() +
      aes(x = x, next_x = next_x, node = node, next_node = next_node,
          fill = node) +
      geom_sankey() +
      scale_fill_manual(values = class_colors) +
      theme_void() +
      theme(legend.position = 'none')
ggsave(filename = 'reassignment_sankey_class.png', path = '05-results/figure1/raw_R_plots/', width = 4, height = 8, dpi = 900)

# subclass sankey
meta |> 
      make_long(subclass_name.with, subclass_name.without) |> 
ggplot() +
      aes(x = x, next_x = next_x, node = node, next_node = next_node,
          fill = node) +
      geom_sankey() +
      scale_fill_manual(values = subclass_colors) +
      theme_void() +
      theme(legend.position = 'none')
ggsave(filename = 'reassignment_sankey_subclass.png', path = '05-results/figure1/raw_R_plots/', width = 4, height = 8, dpi = 900)

# supertype sankey
meta |> 
      filter(subclass_name.with == '016 CA1-ProS Glut') |> 
      make_long(supertype_name.with, supertype_name.without) |> 
ggplot() +
      aes(x = x, next_x = next_x, node = node, next_node = next_node,
          fill = node) +
      geom_sankey() +
      scale_fill_manual(values = supertype_colors) +
      theme_void() +
      theme(legend.position = 'none')
ggsave(filename = 'reassignment_sankey_supertype.png', path = '05-results/figure1/raw_R_plots/', width = 4, height = 8, dpi = 900)

# cluster sankey
meta |> 
      filter(subclass_name.with == '016 CA1-ProS Glut') |> 
      make_long(cluster_name.with, cluster_name.without) |> 
ggplot() +
      aes(x = x, next_x = next_x, node = node, next_node = next_node,
          fill = node) +
      geom_sankey() +
      scale_fill_manual(values = cluster_colors) +
      theme_void() +
      theme(legend.position = 'none')
ggsave(filename = 'reassignment_sankey_cluster.png', path = '05-results/figure1/raw_R_plots/', width = 4, height = 8, dpi = 900)
