# load libs and data ----
print("Loading libraries and data...")
library(Seurat)
library(dplyr)
library(readr)
library(readxl)
library(DESeq2)
library(patchwork)
library(glue)
library(tidyr)
library(ggplot2)
library(plotly)

source('03-scripts/R/seq_functions.R')
nuclei <- readRDS('04-analysis/Seurats/Dec2024/seurat.Rds')
activity_colors <- LoadActivityColors("Dec2024")
subclass_colors <- LoadAllenColors("subclass")


# establish subclasses to use ----
print('Subsetting Seurat object...')
subclass_list <- LoadSubclassesToUse(nuclei, ascertainment = 'custom')
subclass_sets <- c(
  excitatory = list(subclass_list[1:7]),
  inhibitory = list(subclass_list[9:16]),
  glia = list(subclass_list[c(8,17:20)])
)
seurat_subsets <- c(
  excitatory = subset(nuclei, subclass_name %in% subclass_sets[[1]]),
  inhibitory = subset(nuclei, subclass_name %in% subclass_sets[[2]]),
  glia       = subset(nuclei, subclass_name %in% subclass_sets[[3]])
)

# loop thru sets of subclasses
plots <- list()
df_exp <- list()
for (i in seq_along(subclass_sets)) {
  set_name <- names(subclass_sets[i])
  subclass_list <- subclass_sets[[i]]
    
  # get list of subclasses to use
  nuclei_subclass <- seurat_subsets[[i]]
  Idents(nuclei_subclass) <- nuclei_subclass$subclass_name
  
  # read in DEGs ----
  print("Getting DEGs for all subclasses and contrasts...")
  contrast_list <- c("EE30m_vs_SE", "EE6h_vs_SE")
  dir_deg <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
  deg_files <- list.files(dir_deg, full.names = TRUE)
  
  # loop thru CSVs and load in DEGs
  df_all_DEGs <- data.frame()  # init empty df
  for (file in deg_files) {
    fname <- basename(file)
    fname <- str_remove(fname, "\\.csv$")
    # Extract the contrast (last two underscore-separated parts in filename)
    contrast <- str_extract(fname, "[^_]+_[^_]+_[^_]+$")
    # Extract the subclass (everything before the contrast in filename)
    subclass <- str_remove(fname, paste0("_", contrast, "$"))
    subclass <- str_replace_all(subclass, "_", " ")
    subclass <- str_sub(subclass, 1, -2)
    
    # skip excluded subclasses or contrasts
    if (!(contrast %in% contrast_list)) {
      next
    } else if (!(subclass %in% subclass_list)) {
      next
    }
    
    # add gene list to df
    deg <- read_csv(file) |> 
      arrange(desc(log2FoldChange)) |> 
      filter(classification != 'no_change')
    # filter(padj < 0.05)
    
    deg$subclass <- subclass
    deg$contrast <- contrast
    
    df_all_DEGs <- rbind(df_all_DEGs, deg)
  }
  
  
  # df_distinct_DEGs ----
  df_distinct_DEGs <- df_all_DEGs |> 
    filter(subclass %in% subclass_list) |> 
    filter(abs(log2FoldChange) >= 0.585) |> 
    mutate(activity_condition = case_when(
      contrast == "EE30m_vs_SE" ~ "EE30m",
      contrast == "EE6h_vs_SE" ~ "EE6h",
      contrast == "KA30m_vs_SE" ~ "KA30m",
      contrast == "KA6h_vs_SE" ~ "KA6h",
      TRUE ~ NA_character_
    )) |>
    filter(!str_detect(contrast, "KA")) |>
    group_by(gene, subclass, contrast) |> 
    summarize(log2FoldChange, activity_condition, .groups = 'drop') |> 
    distinct(gene, .keep_all = TRUE) |> 
    mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) |> 
    select(-subclass) |>  # remove for future joins
    print()
  
  
  # normalize data DESeq2 ----
  # print("Normalizing expression data...")
  # pseudobulk_counts <- AggregateExpression(
  #   nuclei_subclass,
  #   assays = 'RNA',
  #   features = df_distinct_DEGs$gene,
  #   group.by = c("subclass_name", "activity_condition")
  # )$RNA |>
  #   as.matrix()
  # 
  # coldata <- data.frame(
  #   group = colnames(pseudobulk_counts),
  #   stringsAsFactors = FALSE
  # )
  # coldata <- coldata |>
  #   separate(group, into = c("subclass", "activity_condition"), sep = "_", remove = F)
  # rownames(coldata) <- coldata$group
  # 
  # dds <- DESeqDataSetFromMatrix(
  #   countData = pseudobulk_counts,
  #   colData = coldata,
  #   design = ~ subclass + activity_condition
  # )
  # 
  # # Set "SE" as the reference level for activity_condition:
  # dds$activity_condition <- relevel(dds$activity_condition, ref = "SE")
  # dds <- DESeq(dds)
  # 
  # # log-transform normalized counts
  # norm_counts <- counts(dds, normalized = TRUE)
  # norm_log2 <- log2(norm_counts + 1)  # Adding a pseudocount to avoid log(0)
  # df_norm <- as.data.frame(norm_log2) |>
  #   rownames_to_column(var = "gene") |>
  #   pivot_longer(cols = -gene, names_to = "group", values_to = "avg_expression_log2")
  # 
  # # re-extract variables
  # df_norm <- df_norm |>
  #   separate(group, into = c("subclass", "activity_condition"), sep = "_") |>
  #   mutate(subclass = str_sub(subclass, 2, -1))
  # 
  # # calculate log2fc_from_SE for each gene x subclass combo
  # df_expression.DESeq <- df_norm |>
  #   filter(!str_detect(activity_condition, 'KA')) |>
  #   group_by(gene, subclass) |>
  #   mutate(log2fc_from_SE = avg_expression_log2 - avg_expression_log2[activity_condition == "SE"]) |>
  #   left_join(df_distinct_DEGs, by = 'gene') |>
  #   group_by(gene, subclass) |>
  #   mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
  #   mutate(activity_condition = factor(activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))) |>
  #   mutate(subclass = factor(subclass, levels = subclass_list)) |>
  #   relocate(gene, subclass, activity_condition, avg_expression_log2, log2fc_from_SE) |>
  #   ungroup()
  # 
  # df_expression <- df_expression.DESeq
  
  
  # normalize data Seurat ----
  comparisons <- list(
      "EE30m_vs_SE" = c("EE30m", "SE"),
      "EE6h_vs_SE" = c("EE6h", "SE")
    )
  
  # Function to calculate log2FC for a given subclass and comparison
  calculate_fc <- function(subclass, comparison_name, ident1, ident2) {
    nuclei_subclass |>
      FoldChange(group.by = 'activity_condition',
                 ident.1 = ident1,
                 ident.2 = ident2,
                 subset.ident = subclass,
                 features = df_distinct_DEGs$gene,
                 fc.name = 'log2FoldChange',
                 base = 2) |>
      rownames_to_column('gene') |>
      mutate(subclass = subclass,
             contrast = comparison_name) |>
      select(gene, log2FoldChange, subclass, contrast, pct.1, pct.2) # Keep necessary columns
  }
  
  # Iterate over all subclasses and comparisons, storing results in a single tibble
  df_subclass_contrasts <- expand_grid(subclass = subclass_list, comparison = names(comparisons)) |> 
    mutate(ident1 = str_extract(comparison, '^[^_]+'),  # str before first underscore
           ident2 = str_extract(comparison, '[^_]+$'))  # str after 2nd underscore
  
  df_raw <- df_subclass_contrasts |>
    pmap_dfr(calculate_fc) |> 
    tibble()
  
  # calculate shrunken log2FC by scaling down log2FC by the percent of cells the gene is expressed in
  df_fc <- df_raw |> 
    mutate(pct.mean = rowMeans(cbind(pct.1, pct.2))) |>
    mutate(log2FoldChange.shrink = log2FoldChange * sqrt(sqrt(pct.mean)), .after = log2FoldChange) |> 
    relocate(gene, subclass, contrast, pct.mean, log2FoldChange, log2FoldChange.shrink) |> 
    arrange(gene) |> 
    mutate(IEG = gene %in% LoadGeneList()) |> 
    print()
  
  plot_ly(
    df_fc |> filter(gene %in% LoadGeneList()),
    x = ~pct.mean, y = ~log2FoldChange - log2FoldChange.shrink, color = ~IEG, 
    text = ~paste(gene, subclass, contrast, log2FoldChange, log2FoldChange.shrink, sep = '<br>'),
    type = 'scatter',
    mode = 'markers',
    marker = list(size = 10, color = 'black', opacity = 0.5, symbol = 'circle-open')
  )
    
  df_fc |> filter(gene %in% LoadGeneList()) |> 
  ggplot() +
    aes(x = pct.mean, y = log2FoldChange - log2FoldChange.shrink) +
    geom_point(shape = 21, size = 2, alpha = 1)
  
  
  df_expression.seurat <- df_fc |>
    separate(contrast, into = c('group1', 'group2'), sep = '_vs_') |>
    select(-c(group1, group2)) |>
    left_join(df_distinct_DEGs, by = 'gene', suffix = c('', '.deseq')) |>
    mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
    relocate(log2FoldChange.deseq, .after = log2FoldChange.shrink) |>
    mutate(activity_condition = factor(activity_condition, levels = c('EE30m', 'EE6h'))) |>
    mutate(subclass = factor(subclass, levels = subclass_list)) |>
    mutate(gene = factor(gene, levels = df_distinct_DEGs$gene)) |>
    print()
  
  
  # read gene classes ----
  df_gene_classes <- read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/0_DEG_classifications.csv")
  df_expression <- df_fc |> 
    left_join(df_gene_classes |> select(gene, classification), by = 'gene') |> 
    left_join(df_distinct_DEGs |> select(gene, direction)) |> 
    filter(.by = gene, max(abs(log2FoldChange.shrink)) > 0.585) |> 
    separate(contrast, into = c('group1', 'group2'), sep = '_vs_') |> 
    mutate(activity_condition = factor(group1, levels = c('EE30m', 'EE6h'))) |>
    mutate(subclass = factor(subclass, levels = subclass_list)) |>
    mutate(gene = factor(gene, levels = df_distinct_DEGs$gene)) |>
    print()
  
  df_exp <- c(df_exp, list(df_expression))
  
  # tau() ----
  tau <- function(subclass_expression_vector, direction) {
  # This function takes a vector of expression values for a single gene 
  # across different subclasses and calculates the tau specificity score 
  # for that gene in every subclass. The tau score originates from
  # Yanai et al. (2005). Bioinformatics:
  # https://doi.org/10.1093/bioinformatics/bti042
  # The tau for a gene in subclass i out of n subclasses is defined as:
  #      tau = sum(1 - log2FC_i / max-min(log2FC)) / (n - 1)
  # 
  # Args:
  #  subclass_expression_vector: A vector of log2FC expression values for a single gene
  #  across different subclasses.
  # 
  # Returns:
  #  tau: A numeric value representing the tau specificity score for the gene across subclasses.
    
    # calculate tau
    if (direction == 'up') {
      # subclass_expression_vector <- subclass_expression_vector - min(subclass_expression_vector)
    } else {  # simply reverse the signs
      subclass_expression_vector <- -subclass_expression_vector - min(subclass_expression_vector)
    }
    
    tau <- sum(1 - subclass_expression_vector / max(subclass_expression_vector)) / 
      (length(subclass_expression_vector) - 1)
    return(tau)
  }
  
  
  # find and plot high-tau genes
  df_looping = data.frame(classification = c('ERG'), activity_condition = c('EE30m'))
  
  df_tau <- df_expression |> 
    filter(direction == 'up') |> 
    filter(classification == 'ERG') |>
    filter(activity_condition == 'EE30m') |>
    mutate(tau = tau(log2FoldChange.shrink, direction[1]), 
           .by = 'gene', 
           .after = 'gene') |> 
    arrange(desc(tau)) |> 
    group_by(gene) |> 
    mutate(tau_subclass = case_when(
      direction == 'up'   ~ subclass[which.max(log2FoldChange.shrink)],
      direction == 'down' ~ subclass[which.min(log2FoldChange.shrink)]),
      .after = tau
    ) |> 
    print()
  
  ##### test debug code
  # tmp <- df_expression |> 
  #   filter(direction == 'up') |> 
  #   filter(classification == 'ERG') |>
  #   filter(activity_condition == 'EE30m') |> 
  #   filter(gene == 'Siah3') |>
  #   pull(log2FoldChange)
  # 
  # df_expression |> 
  #   filter(direction == 'up') |> 
  #   filter(classification == 'ERG') |>
  #   filter(activity_condition == 'EE30m') |> 
  #   group_by(gene) |> 
  #   mutate(tau = tau(log2FoldChange.shrink, 'up'), .after = 'gene') |> 
  #   arrange(desc(tau))
  ##### test debug code
  
  # find the highest-tau gene for each subclass
  df_high_tau <- df_tau |> 
    group_by(gene) |> 
    slice_head() |> 
    ungroup() |> 
    slice_max(tau, n = 1, by = tau_subclass) |> 
    mutate(tau_subclass = factor(tau_subclass, levels = subclass_list)) |>
    arrange(tau_subclass) |> 
    print()
  
  p1 <- df_expression |> 
    filter(gene %in% df_high_tau$gene) |>
    filter(classification == df_looping$classification[k]) |>
    filter(activity_condition == df_looping$activity_condition[k]) |> 
    group_by(subclass) |> 
    mutate(log2FC_z = (log2FoldChange.shrink - mean(log2FoldChange.shrink)) / sd(log2FoldChange.shrink)) |> 
    mutate(gene = factor(gene, levels = df_high_tau$gene)) |> 
  ggplot() +
    aes(x = gene, y = subclass, fill = log2FC_z) +
    geom_tile() +
    scale_fill_viridis_c() +
    labs(title = glue("{df_looping$classification[k]} for {names(subclass_sets[i])}")) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 20),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_text(size = 20, angle = 45, hjust = 1),
      axis.text.y = element_text(size = 15),
      # legend.position = 'none'
    )
  print(p1)
  
  p2 <- df_expression |> 
    filter(gene %in% df_high_tau$gene) |>
    filter(classification == df_looping$classification[k]) |>
    filter(activity_condition == df_looping$activity_condition[k]) |> 
    group_by(gene) |> 
    mutate(log2FC_z = (log2FoldChange - mean(log2FoldChange)) / sd(log2FoldChange)) |> 
    mutate(gene = factor(gene, levels = df_high_tau$gene)) |> 
    left_join(df_high_tau |> select(gene, tau_subclass), by = 'gene') |>
  ggplot() +
    aes(x = subclass, y = log2FoldChange, color = tau_subclass, group = gene) +
    geom_hline(yintercept = 0) +
    geom_line(linewidth = 2) +
    scale_color_manual(values = subclass_colors) +
    facet_wrap(~gene, ncol = 1, scales = 'free_y') +
    theme_minimal() +
    labs(title = glue("{df_looping$classification[k]} for {names(subclass_sets[i])}")) +
    theme(
      plot.title = element_text(size = 20),
      strip.text.x = element_text(size = 15),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size = 12, angle = 90, hjust = 1),
      axis.text.y = element_text(size = 15),
      legend.position = 'none',
    )
  print(p2)
  
  plots <- c(plots, list(p1, p2))
  # } # end ERG/LRG loop
} # end subclass set loop




# VlnPlots -----
subclass_list <- subclass_sets[[i]]
genes <- c(
'Gm12940', 'Pde10a', 'Arid5b', 'Bop1'  # no shrinkage
)
plotting_gene <- 'Rn7sk'
VlnPlot(seurat_subsets$glia, group.by = 'subclass_name', split.by = 'activity_condition', plotting_gene) +
  scale_fill_manual(values = activity_colors)
