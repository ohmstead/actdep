# This script classifies all DEGs from EE  conditions found across each gigaclass
# (see figure NEWT for gigaclasses). ERG and LRG labels are mutually exclusive.
# 
# Steps for classification:
# 1. Load Seurat object and subset to excitatory, inhibitory, or glia
# 2. Load DEGs for all subclasses in gigaclass and contrasts
# 3. Provisionally classify as ERG/LRG/both based on when it is DE
# 4. Normalize genome-wide expression data using DESeq2
# 5. Calculate log2FoldChange for each gene x subclass combo
# 6. For genes in "both" category, perform final classification:
#     6a: Force to ERG or LRG based on peak fold-change
#     6b: For genes that may be ERG in one subclass and LRG in another, 
#         force to modal classification
# 7. Save classifications to CSV for that gigaclass

source('03-scripts/R/seq_functions.R')

library(Seurat)
library(DESeq2)


# 1. Load gigaclasses ----------------------------------------
print('Load and subset Seurat object...')
seurat_gigaclasses <- LoadDataset("Dec2024", as_gigaclasses = T)
gigaclasses <- LoadSubclassesToUse(as_gigaclasses = T)

threshold <- 0.585


# 2. Load DEGs ----------------------------------------
print("Getting DEGs for all subclasses and contrasts...")
contrast_list <- c("EE30m_vs_SE", "EE6h_vs_SE")
dir_deg <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk"
deg_files <- list.files(dir_deg, full.names = TRUE)

for (gigaclass in names(gigaclasses)) {  # GIGACLASS LOOP
  subclass_list <- gigaclasses[[gigaclass]]
  nuclei_gigaclass <- seurat_gigaclasses[[gigaclass]]
  
  df_all_DEGs <- data.frame()  # init
  for (file in deg_files) {  # DEG LOADING LOOP
    fname <- basename(file)
    fname <- str_remove(fname, "\\.csv$")
    if (str_detect(fname, '^0_')) {next} # skip irrelevant files
    
    # Extract the contrast (last two underscore-separated parts in filename)
    contrast <- str_extract(fname, "[^_]+_[^_]+_[^_]+$")
    
    # Extract the subclass (everything before the contrast in filename)
    subclass <- str_remove(fname, paste0("_", contrast, "$"))
    subclass <- str_replace_all(subclass, "_", " ")
    subclass <- str_sub(subclass, 1, -2)
    
    # skip irrelevant subclasses and contrasts
    if (!(subclass %in% subclass_list)) {next}
    if (!(contrast %in% contrast_list)) {next}
    
    # add gene list to df
    deg <- read_csv(file, show_col_types = F) |> 
      filter(classification != 'no_change')
    
    deg$DE_ascertainment_subclass <- subclass
    deg$subclass <- subclass
    deg$contrast <- contrast
    
    df_all_DEGs <- rbind(df_all_DEGs, deg) |> 
      arrange(padj)
  }
    
  
  # 3. Classify ERG/LRG/both ----------------------------------------
  print("Filtering DEGs")
  df_gene_classifications <- df_all_DEGs |> 
    mutate(is_ERG = contrast == "EE30m_vs_SE",
           is_LRG = contrast == 'EE6h_vs_SE') |> 
    group_by(gene) |>
    summarize(
      is_ERG = any(is_ERG, na.rm = T),
      is_LRG = any(is_LRG, na.rm = T),
    ) |> 
    mutate(classification = case_when(
      is_ERG & is_LRG ~ "both",
      is_ERG ~ "ERG",
      is_LRG ~ "LRG",
      TRUE ~ NA_character_
    ))
  df_gene_classifications
  
  # for each subclass, get expression for each activity_condition
  df_distinct_DEGs <- df_all_DEGs |> 
    filter(subclass %in% subclass_list) |> 
    filter(abs(log2FoldChange.shrink) >= threshold) |> 
    mutate(activity_condition = case_when(
      contrast == "EE30m_vs_SE" ~ "EE30m",
      contrast == "EE6h_vs_SE" ~ "EE6h",
      contrast == "KA30m_vs_SE" ~ "KA30m",
      contrast == "KA6h_vs_SE" ~ "KA6h",
      TRUE ~ NA_character_
    )) |>
    filter(!str_detect(contrast, "KA")) |>
    group_by(gene, subclass, contrast) |> 
    summarize(log2FoldChange.shrink = log2FoldChange.shrink,
              DE_ascertainment_subclass = dplyr::first(DE_ascertainment_subclass),
              .groups = 'drop') |> 
    left_join(df_gene_classifications, by = 'gene') |> 
    mutate(subclass_by_contrast = paste(subclass, contrast, sep = " x ")) |> 
    distinct(gene, .keep_all = TRUE) |> 
    mutate(direction = ifelse(log2FoldChange.shrink > 0, 'up', 'down')) |> 
    select(-subclass) # remove subclass for future joins
  df_distinct_DEGs
  
  # SEURAT LOG2FC ARCHIVAL CODE
  # comparisons <- list(
  #   "EE30m_vs_SE" = c("EE30m", "SE"),
  #   "EE6h_vs_SE" = c("EE6h", "SE")
  # )
  # 
  # # Function to calculate log2FC for a given subclass and comparison
  # calculate_fc <- function(subclass, comparison_name, ident1, ident2) {
  #   nuclei_subclass |>
  #     FoldChange(group.by = 'activity_condition',
  #                ident.1 = ident1,
  #                ident.2 = ident2,
  #                subset.ident = subclass,
  #                features = df_distinct_DEGs$gene,
  #                fc.name = 'log2FoldChange',
  #                base = 2) |>
  #     rownames_to_column('gene') |>
  #     mutate(subclass = subclass,
  #            contrast = comparison_name) |>
  #     select(gene, log2FoldChange, subclass, contrast) # Keep necessary columns
  # }
  # 
  # # Iterate over all subclasses and comparisons, storing results in a single tibble
  # df_fc <- expand_grid(subclass = subclass_list, comparison = names(comparisons)) |>
  #   mutate(ident1 = map_chr(comparison, ~ comparisons[[.x]][1]),
  #          ident2 = map_chr(comparison, ~ comparisons[[.x]][2])) |>
  #   pmap_dfr(~ calculate_fc(..1, ..2, ..3, ..4)) |>
  #   tibble()
  # 
  # df_expression.seurat <- df_fc |>
  #   separate(contrast, into = c('group1', 'group2'), sep = '_vs_') |>
  #   mutate(activity_condition = factor(group1, levels = c('EE30m', 'EE6h'))) |>
  #   select(-c(group1, group2)) |>
  #   mutate(subclass = factor(subclass, levels = subclass_list)) |>
  #   mutate(gene = factor(gene, levels = df_distinct_DEGs$gene)) |>
  #   left_join(df_distinct_DEGs, by = 'gene') |>
  #   mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
  #   print()
  
  
  # 4. DESeq2-norm counts ----------------------------------------
  print("Normalizing expression data...")
  
  # get gene list
  mat <- nuclei_gigaclass[['SCT']]@counts
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- c(rownames(mat), df_distinct_DEGs$gene) |> unique()
  
  pseudobulk_counts <- AggregateExpression(
    nuclei_gigaclass,
    assays = 'RNA',
    features = gene_list,
    group.by = c("subclass_name", "activity_condition")
  )$RNA
  
  coldata <- data.frame(
    group = colnames(pseudobulk_counts),
    stringsAsFactors = FALSE
  )
  coldata <- coldata |>
    separate(group, into = c("subclass", "activity_condition"), sep = "_", remove = F)
  rownames(coldata) <- coldata$group
  
  dds <- DESeqDataSetFromMatrix(
    countData = pseudobulk_counts,
    colData = coldata,
    design = ~ subclass + activity_condition
  )
  dds <- estimateSizeFactors(dds)
  
  
  # 5. Find pseudobulk FC ----------------------------------------
  norm_counts <- counts(dds, normalized = TRUE)
  norm_log2 <- log2(norm_counts + 1)  # Adding a pseudocount to avoid log(0)
  df_norm <- as_tibble(norm_log2, rownames = 'gene') |>
    filter(gene %in% df_distinct_DEGs$gene) |>
    pivot_longer(cols = -gene, names_to = "group", values_to = "avg_expression_log2") |> 
    separate(group, into = c("subclass", "activity_condition"), sep = "_") |>
    mutate(subclass = str_sub(subclass, 2, -1))
  df_norm
  
  # calculate log2FoldChange for each gene x subclass combo
  df_expression.DESeq <- df_norm |>
    filter(!str_detect(activity_condition, 'KA')) |>
    group_by(gene, subclass) |>
    mutate(log2FC_calculated = avg_expression_log2 - avg_expression_log2[activity_condition == "SE"]) |> 
    left_join(df_distinct_DEGs, by = 'gene') |> 
    group_by(gene, subclass) |>
    mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
    mutate(activity_condition = factor(activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))) |>
    relocate(gene, classification, subclass, activity_condition, log2FC_calculated) |>
    ungroup()
  df_expression.DESeq
  
  
  # 6a: Force ERG/LRG ----------------------------------------
  # note: this is done based on peak FC timepoint (30m or 6h)
  df_expression <- df_expression.DESeq
  expression_range <- df_expression |>
    group_by(gene, subclass) |>
    summarise(
      max_expr_condition = activity_condition[which.max(log2FC_calculated)],
      min_expr_condition = activity_condition[which.min(log2FC_calculated)]
    )
  
  # re-classify as ERG or LRG
  df_expression <- df_expression |>
    left_join(expression_range, by = c("gene", "subclass")) |>
    mutate(classification = if_else(
      classification == "both",
      if_else(
        direction == "up",
        if_else(max_expr_condition == "EE30m", "ERG", "LRG"), # if upregulated
        if_else(min_expr_condition == "EE30m", "ERG", "LRG"), # if downregulated
      ),
      classification  # For genes not labeled "both", keep original classification.
    ),
  )
  df_expression
  
  
  # 6b: Force modal classification ----------------------------------------
  df_expression <- df_expression |>
    group_by(gene) |>
    mutate(modal_classification = names(sort(table(classification), decreasing = TRUE))[1]) |>
    ungroup() |>
    mutate(classification = modal_classification) |>
    select(-modal_classification)
  df_expression
  
  
  # 7: Save classifications ----------------------------------------
  save_dir <- "04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/"
  path_df_classification <- glue("{save_dir}/0_DEG_classifications_{gigaclass}.csv")
  path_df_expression     <- glue("{save_dir}/0_df_expression_{gigaclass}.csv")
  
  df_expression <- df_expression |> 
    relocate(gene, classification, DE_ascertainment_subclass, direction, 
             subclass, activity_condition, log2FC_calculated, log2FoldChange.shrink)
  
  df_classification <- df_expression |> 
    group_by(gene) |> 
    filter(subclass == DE_ascertainment_subclass) |>
    filter(classification == 'ERG' & activity_condition == 'EE30m' |
           classification == 'LRG' & activity_condition == 'EE6h')
  df_classification
  
  # write csv
  df_classification |> write_csv(path_df_classification)
  df_expression     |> write_csv(path_df_expression)
  
  # print summary
  print(glue('DEGs for {gigaclass}:'))
  print(glue('  ERG: {sum(df_classification$classification == "ERG")}'))
  print(glue('  LRG: {sum(df_classification$classification == "LRG")}'))
} # end gigaclass loop
