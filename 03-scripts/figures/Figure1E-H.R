## ---- Fig1E-H
source("03-scripts/R/seq_functions.R")

library(ggsankey)
library(scCustomize)
library(reticulate)

reticulate::use_condaenv("sc")
plot_save_dir <- "05-results/Figure1/raw_R_plots"

nuclei_with <- LoadDataset("Dec2024")

class_colors <- LoadAllenColors('class')
subclass_colors <- LoadAllenColors('subclass')
supertype_colors <- LoadAllenColors('supertype')
cluster_colors <- LoadAllenColors('cluster')

# 1. Run Dec2024 with ARGs through MapMyCells.
      # already done

# 2. Remove ARGs from Seurat object.
gene_list_to_rm <- unique(c(LoadGeneList('tyssowski'), LoadGeneList('IEG')))
features_with <- Features(nuclei_with)
features_without <- features_with[!features_with %in% gene_list_to_rm]

# use DietSeurat to get just the RNA assay, which we will then subset
DefaultAssay(nuclei_with) <- 'RNA'
nuclei_without <- DietSeurat(nuclei_with, 
                             assays = 'RNA', 
                             layers = 'counts', 
                             features = features_without)

# Run Dec2024 without ARGs through MapMyCells ---- 
working_save_dir <- "04-analysis/cluster_reassignment"
# as.anndata(x = nuclei_without, file_path = working_save_dir, file_name = "Dec2024_withoutARGs", main_layer = "counts", other_layers = NULL)


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
      filter(subclass_name.with %in% subclasses_to_use) |> 
      select(-starts_with('emptyDrops_'))


# sankeys ----------------------------------------
plots <- list()

# class
p1 <- meta |> 
  make_long(class_name.with, class_name.without) |> 
  mutate(x = fct_recode(
    x,
    "+ARG" = "class_name.with",
    "-ARG" = "class_name.without",
  )) |>
  mutate(next_x = fct_recode(
    next_x,
    "-ARG" = "class_name.without"
  )) |> 
ggplot() +
  aes(x = x, next_x = next_x, node = node, next_node = next_node,
      fill = node) +
  geom_sankey() +
  scale_fill_manual(values = class_colors) +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank()
    ) +
  labs(
    title = 'Class',
    x = '',
    y = ''
  )
# print(p1)
plots <- c(plots, list(p1))

if (exists("SAVE_PLOTS") & SAVE_PLOTS==TRUE){
    ggsave(plot = p,
          filename = 'reassignment_sankey_class.png', 
          path = plot_save_dir, 
          width = 4, height = 8, dpi = 900)
    svgsave(plot = p + LoadBarebonesTheme(),
          filename = 'reassignment__sankey_class.svg', 
          path = plot_save_dir, 
          width = 4, height = 8)
}

# subclass
p2 <- meta |> 
  make_long(subclass_name.with, subclass_name.without) |> 
  mutate(x = fct_recode(
    x,
    "+ARG" = "subclass_name.with",
    "-ARG" = "subclass_name.without",
  )) |>
  mutate(next_x = fct_recode(
    next_x,
    "-ARG" = "subclass_name.without"
  )) |> 
ggplot() +
  aes(x = x, next_x = next_x, node = node, next_node = next_node,
      fill = node) +
  geom_sankey() +
  scale_fill_manual(values = subclass_colors) +
  theme(
    legend.position = 'none',
    axis.text.y = element_blank()
    ) +
  labs(
    title = 'Subclass',
    x = '',
    y = ''
  )
# print(p2)
plots <- c(plots, list(p2))

if (exists("SAVE_PLOTS") & SAVE_PLOTS==TRUE){
    ggsave(plot = p,
          filename = 'reassignment_sankey_subclass.png', 
          path = plot_save_dir, 
          width = 4, height = 8, dpi = 900)
    svgsave(plot = p + LoadBarebonesTheme(),
          filename = 'reassignment__sankey_subclass.svg', 
          path = plot_save_dir, 
          width = 4, height = 8)
}

subclass_list <- list(
  CA1 = '016 CA1-ProS Glut',
  DG  = '037 DG Glut'
)

for (celltype in names(subclass_list)) {
  # supertype sankey
  p <- meta |> 
    filter(subclass_name.with == subclass_list[[celltype]]) |>
    make_long(supertype_name.with, supertype_name.without) |> 
    mutate(x = fct_recode(
      x,
      "+ARG" = "supertype_name.with",
      "-ARG" = "supertype_name.without",
    )) |>
    mutate(next_x = fct_recode(
      next_x,
      "-ARG" = "supertype_name.without"
    )) |> 
  ggplot() +
    aes(x = x, next_x = next_x, node = node, next_node = next_node,
        fill = node) +
    geom_sankey() +
    scale_fill_manual(values = supertype_colors) +
    theme(
      legend.position = 'none',
      axis.text.y = element_blank()
      ) +
    labs(
      title = glue('{celltype} Supertype'),
      x = '',
      y = ''
    )
  plots <- c(plots, list(p))
  
  if (exists("SAVE_PLOTS") & SAVE_PLOTS==TRUE){
      ggsave(plot = p,
            filename = glue('reassignment_sankey_supertype_{celltype}.png'), 
            path = plot_save_dir, 
            width = 4, height = 8, dpi = 900)
      svgsave(plot = p + LoadBarebonesTheme(),
            filename = glue('reassignment__sankey_supertype_{celltype}.svg'), 
            path = plot_save_dir, 
          width = 4, height = 8)
  }
  
  # cluster
  p <- meta |> 
    filter(subclass_name.with == subclass_list[[celltype]]) |>
    make_long(cluster_name.with, cluster_name.without) |> 
    mutate(x = fct_recode(
      x,
      "+ARG" = "cluster_name.with",
      "-ARG" = "cluster_name.without",
    )) |>
    mutate(next_x = fct_recode(
      next_x,
      "-ARG" = "cluster_name.without"
    )) |> 
  ggplot() +
    aes(x = x, next_x = next_x, node = node, next_node = next_node,
        fill = node) +
    geom_sankey() +
    scale_fill_manual(values = cluster_colors) +
    theme(
      legend.position = 'none',
      axis.text.y = element_blank()
      ) +
    labs(
      title = glue('{celltype} Cluster'),
      x = '',
      y = ''
    )
  # print(p)
  plots <- c(plots, list(p))

  if (exists("SAVE_PLOTS") & SAVE_PLOTS==TRUE){
      ggsave(plot = p,
            filename = glue('reassignment_sankey_cluster_{celltype}.png'),
            path = plot_save_dir, 
            width = 4, height = 8, dpi = 900)
      svgsave(plot = p + LoadBarebonesTheme(),
            filename = glue('reassignment__sankey_cluster_{celltype}.svg'), 
            path = plot_save_dir, 
            width = 4, height = 8)
  }
}

# put all plots in list using plot_layout
combined <- Reduce(`+`, plots)
combined + plot_layout(nrow=1)
## ----