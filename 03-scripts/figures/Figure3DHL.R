## ----Fig3DHL
source('03-scripts/R/seq_functions.R')

activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")
save_path <- "05-results/Figure3/raw_R_plots"


# establish subclasses to use ----------------------------------------
print('Subsetting Seurat object...')
gigaclasses <- LoadSubclassesToUse(nuclei, as_gigaclasses = T)
if(!exists('seurat_gigaclass')) {seurat_gigaclass <- LoadDataset("Dec2024", as_gigaclasses = T)}

ieg_symbols <- c(LoadGeneList('IEG'), 'Rn7sk', 'Midn', 
                 'Etv1', 'Rcan2', 'Hs3st2',
                 'Tac1', 'Socs2', 'Daam2', 'Akap5', 'Frmd6', 'Crhbp', 'Gm49673',
                 'Map3k19', 'Usp53', 'Slco1c1', 'Kcnn2')

my_theme <- theme_bw() +
  theme(
    legend.position = 'none',
    plot.title = element_blank(),
    axis.title = element_blank(),
    axis.text.x = element_blank(),
  )
plots <- list()

# excitatory plots ----------------------------------------
genes_to_plot <- c('Hs3st2', 'Rcan2')
df_expression <- read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/0_df_expression_excitatory.csv")
for (gene_to_plot in genes_to_plot) {
  # make ggplots
  p1 <- df_expression |> 
    filter(gene == gene_to_plot, activity_condition == 'EE30m') |> 
    mutate(subclass = factor(subclass, levels = c(gigaclasses$excitatory))) |> 
  ggplot() +
    aes(x = subclass, y = log2FC_calculated, group = gene) +
    geom_line(linewidth = 2) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(title = gene_to_plot) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=30, hjust=1))
  # print(p1)
  
  # make VlnPlot
  seurat_gigaclass$excitatory$subclass_name <- factor(seurat_gigaclass$excitatory$subclass_name, 
                                                  levels = c(gigaclasses$excitatory))
  p2 <- seurat_gigaclass$excitatory |> 
    subset(activity_condition %in% c('SE','EE30m','EE6h')) |> 
    VlnPlot(gene_to_plot, group.by = 'subclass_name', split.by = 'activity_condition') +
    scale_fill_manual(values = activity_colors) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle=30, hjust=1),
      axis.title = element_blank()
    )
  # print(p2)
  plots <- c(plots, list(p2))
  
  # png
  if (exists('SAVE_PLOTS') & SAVE_PLOTS == T) {
    ggsave(plot = p1,
           filename = glue("line_excitatory_{gene_to_plot}.png"),
           path = save_path,
           width = 4, height = 2)
    ggsave(plot = p2,
           filename = glue("VlnPlot_excitatory_{gene_to_plot}.png"),
           path = save_path,
           width = 4, height = 2)
    
    # svg
    ggsave(plot = p1 + LoadBarebonesTheme(ticks = 'y'),
           filename = glue("line__excitatory_{gene_to_plot}.svg"),
           path = save_path,
           width = 4, height = 2)
    ggsave(plot = p2 + LoadBarebonesTheme(ticks = 'y'),
           filename = glue("VlnPlot__excitatory_{gene_to_plot}.svg"),
           path = save_path,
           width = 4, height = 2)
  }
}


# inhibitory plots ----------------------------------------
genes_to_plot <- c('Tac1', 'Crhbp')
df_expression <- read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/0_df_expression_inhibitory.csv")
for (gene_to_plot in genes_to_plot) {
  # make ggplots
  p1 <- df_expression |> 
    filter(gene == gene_to_plot, activity_condition == 'EE30m') |> 
    mutate(subclass = factor(subclass, levels = c(gigaclasses$inhibitory))) |> 
    ggplot() +
    aes(x = subclass, y = log2FC_calculated, group = gene) +
    geom_line(linewidth = 2) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(title = gene_to_plot) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=30, hjust=1))
  # print(p1)
  
  # make VlnPlot
  seurat_gigaclass$inhibitory$subclass_name <- factor(seurat_gigaclass$inhibitory$subclass_name, 
                                                      levels = c(gigaclasses$inhibitory))
  p2 <- seurat_gigaclass$inhibitory |> 
    subset(activity_condition %in% c('SE','EE30m','EE6h')) |> 
    VlnPlot(gene_to_plot, group.by = 'subclass_name', split.by = 'activity_condition') +
    scale_fill_manual(values = activity_colors) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle=30, hjust=1),
      axis.title = element_blank()
    )
  # print(p2)
  plots <- c(plots, list(p2))
  
  # png
  if (exists('SAVE_PLOTS') & SAVE_PLOTS == T) {
    ggsave(plot = p1,
           filename = glue("line_inhibitory_{gene_to_plot}.png"),
           path = save_path,
           width = 4, height = 2)
    ggsave(plot = p2,
           filename = glue("VlnPlot_inhibitory_{gene_to_plot}.png"),
           path = save_path,
           width = 4, height = 2)
    
    # svg
    ggsave(plot = p1 + LoadBarebonesTheme(ticks = 'y'),
           filename = glue("line__inhibitory_{gene_to_plot}.svg"),
           path = save_path,
           width = 4, height = 2)
    ggsave(plot = p2 + LoadBarebonesTheme(ticks = 'y'),
           filename = glue("VlnPlot__inhibitory_{gene_to_plot}.svg"),
           path = save_path,
           width = 4, height = 2)
    }
}


# glia plots ----------------------------------------
genes_to_plot <- c('Map3k19', 'Slco1c1')
df_expression <- read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/0_df_expression_glia.csv")
for (gene_to_plot in genes_to_plot) {
  # make ggplots
  p1 <- df_expression |> 
    filter(gene == gene_to_plot, activity_condition == 'EE30m') |> 
    mutate(subclass = factor(subclass, levels = c(gigaclasses$glia))) |> 
    ggplot() +
    aes(x = subclass, y = log2FC_calculated, group = gene) +
    geom_line(linewidth = 2) +
    geom_hline(yintercept = 0, linetype = 'dashed') +
    labs(title = gene_to_plot) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=30, hjust=1))
  # print(p1)
  
  # make VlnPlot
  seurat_gigaclass$glia$subclass_name <- factor(seurat_gigaclass$glia$subclass_name, 
                                                      levels = c(gigaclasses$glia))
  p2 <- seurat_gigaclass$glia |> 
    subset(activity_condition %in% c('SE','EE30m','EE6h')) |> 
    VlnPlot(gene_to_plot, group.by = 'subclass_name', split.by = 'activity_condition') +
    scale_fill_manual(values = activity_colors) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle=30, hjust=1),
      axis.title = element_blank()
    )
  # print(p2)
  plots <- c(plots, list(p2))
  
  # png
  if (exists('SAVE_PLOTS') & SAVE_PLOTS == T) {
    ggsave(plot = p1,
           filename = glue("line_glia_{gene_to_plot}.png"),
           path = save_path,
           width = 4, height = 2)
    ggsave(plot = p2,
           filename = glue("VlnPlot_glia_{gene_to_plot}.png"),
           path = save_path,
           width = 4, height = 2)
    
    # svg
    ggsave(plot = p1 + LoadBarebonesTheme(ticks = 'y'),
           filename = glue("line__glia_{gene_to_plot}.svg"),
           path = save_path,
           width = 4, height = 2)
    ggsave(plot = p2 + LoadBarebonesTheme(ticks = 'y'),
           filename = glue("VlnPlot__glia_{gene_to_plot}.svg"),
           path = save_path,
           width = 4, height = 2)
  }
}

p <- Reduce('+', plots) + plot_layout(ncol = 2)
print(p)
## ----