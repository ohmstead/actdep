## ----Fig5A
# This script produces simple plots of Clock genes at different ZTs
# as part of figure MOTH.

source('03-scripts/R/seq_functions.R')

library(rstatix)

if (!exists('nuclei')) {nuclei <- LoadDataset('Dec2024')}

activity_colors <- LoadActivityColors('Dec2024')
subclass_colors <- LoadAllenColors('subclass')
subclass_list <- c(
  '016 CA1-ProS Glut',
  '037 DG Glut',
  '319 Astro-TE NN',
  '327 Oligo NN'
)

nuclei_subclass <- nuclei |> 
  subset(activity_condition == 'SE' & subclass_name %in% subclass_list)


calc_sem <- function(data, varname, groupnames) {
  #+++++++++++++++++++++++++
  # Calculatse mean and sd for each group
  #+++++++++++++++++++++++++
  # data : a data frame
  # varname : the name of a column containing the variable
  #           to be summarized
  # groupnames : vector of column names to be used as
  #              grouping variables
  require(plyr)
  
  # Step 1: Summarize expression within each ZT, subclass, and gene group
  summary_func <- function(x) {
    # Count non-missing values
    n <- sum(!is.na(x$expression))
    
    # Calculate mean, standard deviation, and SEM for expression
    mean <- mean(x$expression, na.rm = TRUE)
    expression_sd   <- sd(x$expression, na.rm = TRUE)
    sem  <- expression_sd / sqrt(n)
    
    c(mean = mean,
      sem  = sem)
  }
  
  # Group by ZT, subclass_name, and gene to compute summary statistics
  data_sum <- ddply(data, c("ZT", "subclass_name", "gene"), summary_func)
  
  # Step 2: Normalize the expression values within each gene and subclass across ZTs
  # Get the maximum expression mean per gene in each subclass
  data_sum <- ddply(data_sum, c("gene", "subclass_name"),
                    transform, 
                    max_expression_mean = max(mean, na.rm = TRUE))
  
  # Create normalized columns by dividing by the maximum mean
  data_sum <- transform(data_sum,
                        mean_norm = mean / max_expression_mean,
                        sem_norm  = sem  / max_expression_mean)
  
  # Remove the temporary column
  data_sum$max_expression_mean <- NULL
  
  return(data_sum)
}


# SE ZT insets ----------------------------------------
known_circadian_genes <- LoadGeneList('circadian')
mat_circadian <- GetAssayData(nuclei_subclass)
mat_circadian <- mat_circadian[known_circadian_genes,] |> t()
nuclei_subclass@meta.data <- nuclei_subclass@meta.data |> 
  mutate(ZT.collection = case_when(
    (ZT == 'ZT0' & activity_condition == 'EE6h') ~ 'ZT6',
    (ZT == 'ZT4' & activity_condition == 'EE6h') ~ 'ZT10',
    (ZT == 'ZT12' & activity_condition == 'EE6h') ~ 'ZT18',
    (ZT == 'ZT16' & activity_condition == 'EE6h') ~ 'ZT22',
    TRUE ~ ZT
  )) |> 
  mutate(ZT.collection = factor(ZT.collection, 
                                levels = c('ZT0', 'ZT4', 'ZT6', 'ZT10', 'ZT12', 'ZT16', 'ZT18', 'ZT22')))

# merge circadian gene expression with cell metadata
df_circadian <- as.matrix(mat_circadian) |> 
  as_tibble(rownames = 'barcode') |> 
  right_join(nuclei_subclass@meta.data) |> 
  pivot_longer(cols = all_of(known_circadian_genes), names_to = 'gene', values_to = 'expression') |> 
  select(gene, barcode, ZT, activity_condition, subclass_name, expression) |> 
  calc_sem('mean_expression', c('ZT', 'subclass_name', 'gene')) |> 
  arrange(gene, subclass_name, ZT) |> 
  relocate(gene, subclass_name) |> 
  as_tibble()


# plot clock genes ----------------------------------------
genes_to_plot <- c('Per1', 'Per2', 'Cry2', 'Clock', 'Bmal1')
plots <- list()

for (current_gene in genes_to_plot) {
  p <- df_circadian |> 
    filter(gene == current_gene) |>
  ggplot() +
    aes(x = ZT, y = mean_norm, group = subclass_name, color = subclass_name) +
    geom_line(linewidth = 2) +
    geom_errorbar(aes(ymin = mean_norm-sem_norm, ymax = mean_norm+sem_norm), 
                  linewidth = 2, width = 0.1) +
    geom_point(size = 5) +
    scale_color_manual(values = subclass_colors) +
    labs(title = glue('{current_gene}'),
         y = 'Normalized expression') +
    theme(
      legend.position = 'none',
      axis.text.x = element_text(size = 12, angle = 30),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_blank(),
      plot.title = element_text(size = 12, hjust = 0.5, face = 'italic')
    )
  plots <- append(plots, list(p))
  
  if (SAVE_PLOTS) {
    save_path <- "05-results/Figure5/raw_R_plots/4_subclass_plots"
    ggsave(plot = p,  # png
           filename = glue('inset__{current_gene}.png'),
           path = save_path,
           width = 4, height = 4, units = 'in', dpi = 900, bg = 'white')
    ggsave(plot = p + LoadBarebonesTheme(ticks = 'y'),  # svg
           filename = glue('inset_{current_gene}.svg'),
           path = save_path,
           width = 5, height = 5, units = 'in')
  }
}

p <- Reduce(`+`, plots) + plot_layout(nrow = 1, )
p
## ----