library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

library(Seurat)
library(ComplexHeatmap)

source('03-scripts/R/seq_functions.R')
activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# establish subclasses to use ----------------------------------------
print('Subsetting Seurat object...')
subclass_sets <- LoadSubclassesToUse(nuclei, as_gigaclasses = T)
# seurat_subsets <- LoadDataset("Dec2024", as_gigaclasses = T)

ieg_symbols <- c(LoadGeneList('IEG'), 'Rn7sk', 'Rcan2', 'Tac1', 'Map3k19')

my_theme <- theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
  )

# Rcan2 excitatory plots ----------------------------------------
gene_to_plot <- 'Rcan2'
df_expression <- read_csv(glue("04-analysis/df_expression/df_expression.seurat/df_expression_excitatory.csv"))
p1 <- df_expression |> 
  filter(gene == gene_to_plot, activity_condition == 'EE30m') |> 
  mutate(subclass = factor(subclass, levels = c(subclass_sets$excitatory))) |> 
ggplot() +
  aes(x = subclass, y = log2FoldChange.shrink, group = gene) +
  geom_line(linewidth = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  my_theme
p1
si(400,200)

seurat_subsets$excitatory$subclass_name <- factor(seurat_subsets$excitatory$subclass_name, 
                                                  levels = c(subclass_sets$excitatory))
p2 <- seurat_subsets$excitatory |> 
  subset(activity_condition %in% c('SE','EE30m','EE6h')) |> 
  VlnPlot(gene_to_plot, group.by = 'subclass_name', split.by = 'activity_condition') +
  scale_fill_manual(values = activity_colors) +
  my_theme
p2
si(800,300)

p1 / p2 + plot_layout(heights = c(1,2))
si(800,500)


# Tac1 inhibitory plots ----------------------------------------
gene_to_plot <- 'Tac1'
df_expression <- read_csv(glue("04-analysis/df_expression/df_expression.seurat/df_expression_inhibitory.csv"))
p1 <- df_expression |> 
  filter(gene == gene_to_plot, activity_condition == 'EE30m') |> 
  mutate(subclass = factor(subclass, levels = c(subclass_sets$inhibitory))) |> 
  ggplot() +
  aes(x = subclass, y = log2FoldChange.shrink, group = gene) +
  geom_line(linewidth = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  my_theme
p1
si(400,200)

seurat_subsets$inhibitory$subclass_name <- factor(seurat_subsets$inhibitory$subclass_name, 
                                            levels = c(subclass_sets$inhibitory))
p2 <- seurat_subsets$inhibitory |> 
  subset(activity_condition %in% c('SE','EE30m','EE6h')) |> 
  VlnPlot(gene_to_plot, group.by = 'subclass_name', split.by = 'activity_condition') +
  scale_fill_manual(values = activity_colors) +
  my_theme
p2
si(800,300)

p1 / p2 + plot_layout(heights = c(1,2))
si(800,500)


# Map3k19 glia plots ----------------------------------------
gene_to_plot <- 'Map3k19'
df_expression <- read_csv(glue("04-analysis/df_expression/df_expression.seurat/df_expression_glia.csv"))
p1 <- df_expression |> 
  filter(gene == gene_to_plot, activity_condition == 'EE30m') |> 
  mutate(subclass = factor(subclass, levels = c(subclass_sets$glia))) |> 
  ggplot() +
  aes(x = subclass, y = log2FoldChange.shrink, group = gene) +
  geom_line(linewidth = 2) +
  geom_hline(yintercept = 0, linetype = 'dashed') +
  theme_void() +
  my_theme

p1
si(400,200)

seurat_subsets$glia$subclass_name <- factor(seurat_subsets$glia$subclass_name, 
                                            levels = c(subclass_sets$glia))
p2 <- seurat_subsets$glia |> 
  subset(activity_condition %in% c('SE','EE30m','EE6h')) |> 
VlnPlot(gene_to_plot, group.by = 'subclass_name', split.by = 'activity_condition') +
  scale_fill_manual(values = activity_colors) +
  my_theme

p2
si(800,300)

p1 / p2 + plot_layout(heights = c(1,2))
si(800,500)
