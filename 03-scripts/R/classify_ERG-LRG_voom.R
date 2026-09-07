# This script classifies all DEGs from EE  conditions found across each gigaclass
# (see figure NEWT for gigaclasses). ERG and LRG labels are mutually exclusive.
# 
# Steps for classification:
# 1. Load Seurat object and subset to excitatory, inhibitory, or glia
# 2. Load DEGs for all subclasses in gigaclass and contrasts
# 3. Provisionally classify as ERG/LRG/both based on when it is DE
# 4. Normalize genome-wide expression data using limma-voom (TMM + log2-CPM)
# 5. Calculate log2FoldChange for each gene x subclass combo
# 6. For genes in "both" category, perform final classification:
#     6a: Force to ERG or LRG based on peak fold-change
#     6b: For genes that may be ERG in one subclass and LRG in another, 
#         force to modal classification
# 7. Save classifications to CSV for that gigaclass

source('03-scripts/R/seq_functions.R')

if (!exists("DEG_METHOD"))   DEG_METHOD   <- "voom"
if (!exists("dir_deg"))      dir_deg      <- "04-analysis/reviewer_response/voom"
if (!exists("save_dir"))     save_dir     <- dir_deg
save_dir <- sub("/+$", "", save_dir)  # normalize any trailing slash(es)

library(Seurat)
library(limma)
library(edgeR)


# 1. Load gigaclasses ----------------------------------------
# NOTE: the three gigaclass objects total ~17 GB in memory, which does not fit
# alongside the pseudobulk matrices on an 18 GB machine. Each one is loaded
# inside the gigaclass loop and dropped again before the next iteration.
gigaclasses <- LoadSubclassesToUse(as_gigaclasses = T)
gigaclass_paths <- c(
  excitatory = "04-analysis/Seurats/Dec2024/seurat_excitatory.Rds",
  inhibitory = "04-analysis/Seurats/Dec2024/seurat_inhibitory.Rds",
  glia       = "04-analysis/Seurats/Dec2024/seurat_glia.Rds"
)

threshold <- 0.585


# 2. Load DEGs ----------------------------------------
print("Getting DEGs for all subclasses and contrasts...")
contrast_list <- c("EE30m_vs_SE", "EE6h_vs_SE")
deg_files <- list.files(dir_deg, full.names = TRUE)
deg_files <- deg_files[!str_detect(basename(deg_files), "^\\._")]  # drop AppleDouble sidecar files

for (gigaclass in names(gigaclasses)) {  # GIGACLASS LOOP
  subclass_list <- gigaclasses[[gigaclass]]
  print(glue("Loading {gigaclass} Seurat object..."))
  nuclei_gigaclass <- LoadSeuratRds(gigaclass_paths[[gigaclass]])
  
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
    deg <- read_csv(file, show_col_types = F)
    deg <- filter(deg, classification != 'no_change')
    
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
    filter(abs(log2FoldChange) >= threshold) |> 
    mutate(activity_condition = case_when(
      contrast == "EE30m_vs_SE" ~ "EE30m",
      contrast == "EE6h_vs_SE" ~ "EE6h",
      contrast == "KA30m_vs_SE" ~ "KA30m",
      contrast == "KA6h_vs_SE" ~ "KA6h",
      TRUE ~ NA_character_
    )) |>
    filter(!str_detect(contrast, "KA")) |>
    group_by(gene, subclass, contrast) |> 
    summarize(log2FoldChange = log2FoldChange,
              DE_ascertainment_subclass = dplyr::first(DE_ascertainment_subclass),
              .groups = 'drop') |> 
    left_join(df_gene_classifications, by = 'gene') |> 
    mutate(subclass_by_contrast = paste(subclass, contrast, sep = " x ")) |> 
    distinct(gene, .keep_all = TRUE) |> 
    mutate(direction = ifelse(log2FoldChange > 0, 'up', 'down')) |> 
    select(-subclass) # remove subclass for future joins
  df_distinct_DEGs
  
  
  # 4. voom-normalized counts ----------------------------------------
  # Normalization mirrors the pseudobulk limma-voom DGE models
  # (manuscript_DGE_reviewer_response_02_pseudobulk_voom.R): animal-level
  # pseudobulk counts -> TMM norm factors -> voom log2-CPM.
  print("Normalizing expression data (limma-voom)...")

  # get gene list
  mat <- nuclei_gigaclass[['SCT']]@counts
  mat <- mat[rowMeans(mat > 0) > 0.01, ]
  gene_list <- c(rownames(mat), df_distinct_DEGs$gene) |> unique()
  gene_list <- intersect(gene_list, rownames(nuclei_gigaclass[['RNA']]))

  # pseudobulk at the animal level, the same unit of analysis as the voom models
  pseudobulk_counts <- AggregateExpression(
    nuclei_gigaclass,
    assays = 'RNA',
    features = gene_list,
    group.by = c("subclass_name", "sample")
  )$RNA
  pseudobulk_counts <- as.matrix(pseudobulk_counts)

  # metadata is all that is still needed; release the Seurat object
  meta_gigaclass <- nuclei_gigaclass@meta.data
  rm(nuclei_gigaclass)
  gc()

  coldata <- data.frame(
    group = colnames(pseudobulk_counts),
    stringsAsFactors = FALSE
  )
  coldata <- coldata |>
    separate(group, into = c("subclass", "sample"), sep = "_", remove = F) |>
    mutate(subclass = str_sub(subclass, 2, -1))  # drop the "g" Seurat prepends

  # map each pseudobulk column back to its activity_condition
  sample_meta <- meta_gigaclass |>
    as_tibble() |>
    distinct(sample, activity_condition) |>
    mutate(sample = str_replace_all(sample, "_", "-"))  # AggregateExpression naming
  coldata <- left_join(coldata, sample_meta, by = "sample")
  rownames(coldata) <- coldata$group
  stopifnot(!any(is.na(coldata$activity_condition)))

  # TMM normalization + voom, with each subclass x condition as its own group
  grp <- factor(paste(coldata$subclass, coldata$activity_condition, sep = " x "))
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- make.names(levels(grp))

  y <- DGEList(counts = pseudobulk_counts, group = grp)
  keep <- filterByExpr(y, design = design) | rownames(y) %in% df_distinct_DEGs$gene
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  v <- voom(y, design)


  # 5. Find pseudobulk FC ----------------------------------------
  # voom log2-CPM, averaged across animals within subclass x activity_condition
  df_norm <- as_tibble(v$E, rownames = 'gene') |>
    filter(gene %in% df_distinct_DEGs$gene) |>
    pivot_longer(cols = -gene, names_to = "group", values_to = "expression_log2") |>
    left_join(dplyr::select(coldata, group, subclass, activity_condition),
              by = "group") |>
    group_by(gene, subclass, activity_condition) |>
    summarize(avg_expression_log2 = mean(expression_log2), .groups = 'drop')
  df_norm

  # calculate log2FoldChange for each gene x subclass combo
  df_expression.voom <- df_norm |>
    filter(!str_detect(activity_condition, 'KA')) |>
    group_by(gene, subclass) |>
    mutate(log2FC_calculated = avg_expression_log2 - avg_expression_log2[activity_condition == "SE"]) |>
    left_join(df_distinct_DEGs, by = 'gene') |>
    group_by(gene, subclass) |>
    mutate(subclass_by_activity_condition = paste(subclass, activity_condition, sep = ' x ')) |>
    mutate(activity_condition = factor(activity_condition, levels = c('SE', 'EE30m', 'EE6h', 'KA30m', 'KA6h'))) |>
    relocate(gene, classification, subclass, activity_condition, log2FC_calculated) |>
    ungroup()
  df_expression.voom


  # 6a: Force ERG/LRG ----------------------------------------
  # note: this is done based on peak FC timepoint (30m or 6h)
  df_expression <- df_expression.voom
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
  path_df_classification <- glue("{save_dir}/0_DEG_classifications_{gigaclass}.csv")
  path_df_expression     <- glue("{save_dir}/0_df_expression_{gigaclass}.csv")
  
  df_expression <- df_expression |> 
    relocate(gene, classification, DE_ascertainment_subclass, direction, 
             subclass, activity_condition, log2FC_calculated, log2FoldChange)
  
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

  rm(pseudobulk_counts, y, v, meta_gigaclass)
  gc()
} # end gigaclass loop
