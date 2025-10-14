# A script housing a collection of commonly called functions in the actdep project.

# helper to attach packages only when needed ------------------------------
attach_package_once <- function(pkg) {
  if (!paste0("package:", pkg) %in% search()) {
    library(pkg, character.only = TRUE)
  }
}

# load common libs ----------------------------------------
packages_to_attach <- c(
  "ggplot2",
  "plotly",
  "patchwork",
  "glue",
  "forcats",
  "readr",
  "readxl",
  "stringr",
  "tibble",
  "tidyr",
  "Seurat"
)
invisible(lapply(packages_to_attach, attach_package_once))

# always load at end of tidyverse to mask plyr; plyr fucking sucks
attach_package_once("dplyr")

# This function will get called every time the script is sourced!
set_ggplot_font <- function(font = 'Aptos') {
  theme_set(theme_minimal(base_family = font))
  message(glue("{font} is now the default font for ggplot2."))
}
set_ggplot_font('Helvetica') # set default font


SymbolToEnsembl <- function(gene_symbol_vector) {
  # This function takes gene symbols and converts them to Ensemble IDs.
  # Useful for GO analysis.
  gene_dict <- read_csv("02-data/raw_data/all_genes.csv", show_col_types = F)
  ensembl_ids <- gene_dict |> 
    filter(gene_name %in% gene_symbol_vector) |> 
    pull(gene_id)
  
  return(ensembl_ids)
}


tau <- function(subclass_expression_vector, direction = 'up') {
  # This function takes a vector of expression values for a single gene 
  # across different subclasses and calculates the tau specificity score 
  # for that gene in every subclass. The tau score originates from
  # Yanai et al. (2005). Bioinformatics:
  # https://doi.org/10.1093/bioinformatics/bti042
  # The tau for a gene in subclass i out of n subclasses is defined as:
  #      tau = sum(1 - log2FC_i / max(log2FC)) / (n - 1)
  # 
  # Args:
  #  subclass_expression_vector: A vector of log2FC expression values for a single gene
  #  across different subclasses.
  # 
  # Returns:
  #  tau: A numeric value representing the tau specificity score for the gene across subclasses.
  
  # switch direction if downregulated gene
  if (direction == 'down') {
    subclass_expression_vector <- -subclass_expression_vector
  }
  
  # calculate tau
  tau <- sum(1 - subclass_expression_vector / max(subclass_expression_vector)) / 
    (length(subclass_expression_vector) - 1)
  return(tau)
}


LoadDataset <- function(dataset, as_gigaclasses = FALSE, sublibrary = "combined") {
  attach_package_once("Seurat")
  
  # try gigaclasses first
  if (as_gigaclasses) {
    load_dir <- glue("04-analysis/Seurats/{dataset}/")
    seurat_gigaclasses <- c(
      excitatory = LoadSeuratRds(paste0(load_dir, "seurat_excitatory.Rds")),
      inhibitory = LoadSeuratRds(paste0(load_dir, "seurat_inhibitory.Rds")),
      glia       = LoadSeuratRds(paste0(load_dir, "seurat_glia.Rds"))
    )
    return(seurat_gigaclasses)  # return list of seurats
  }
  
  # try Yao next
  if (dataset %in% c('yao', 'Yao', 'yao2023', 'Yao2023')) {
    # try to load from local file
    if (file.exists("04-analysis/Seurats/Yao2023/seurat.Rds")){
      local_file_path <- "04-analysis/Seurats/Yao2023/seurat.Rds"
    } else {
      data_path <- "/Volumes/jack/Yao2023/"
      seurat_obj <- LoadSeuratRds(paste0(data_path, "yao_seurat.rds"))
    }
  } else {  # try e.g. "Dec2024"
    local_file_path <- glue("04-analysis/Seurats/{dataset}/seurat.Rds")
    ssd_file_path <- glue("/Volumes/jack/seq/analysis_{dataset}/{sublibrary}/0_all-sample/DGE_filtered/seurat.Rds")
    if (file.exists(local_file_path)) {
      seurat_obj <- LoadSeuratRds(local_file_path)
    } else {
      seurat_obj <- LoadSeuratRds(ssd_file_path)
    }
  }
}


LoadAllenColors <- function(clade = 'subclass') {
  allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv", show_col_types = FALSE)
  
  switch(
    clade,
    'class' = return(allen_colors |> distinct(class, class_color) |> deframe()),
    'subclass' = return(allen_colors |> distinct(subclass, subclass_color) |> deframe()),
    'supertype' = return(allen_colors |> distinct(supertype, supertype_color) |> deframe()),
    'cluster' = return(allen_colors |> distinct(cluster, cluster_color) |> deframe())
  )
}


LoadGeneList <- function(list_type = "IEG") {
  # has options:
  #   - IEG
  #   - lncRNA: lncRNA genes from CA1 EE30m vs SE
  #   - tyssowski: tyssowski rapid/delayed PRGs
  #   - DEGs: all DEGs from the expression heatmap in figure 1
  #   - circadian: Clock genes
  attach_package_once("readxl")
  attach_package_once("dplyr")
  
  if (list_type == "IEG" | list_type == 'ieg') {
    gene_list <- c(
      'Arc',
      'Btg2',
      'Dusp1',
      'Dusp5',
      'Egr1',
      'Egr2',
      'Egr3',
      'Egr4',
      'Fos',
      'Fosb',
      'Fosl2',
      'Junb',
      'Npas4',
      'Nr4a1',
      'Nr4a3'
    )
  } else if (list_type == "lncRNA") {
    gene_biotypes <- read_csv("04-analysis/all_gene_stats.csv")
    gene_list <- c(
      read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/016_CA1-ProS_Glut__EE30m_vs_SE.csv") |> 
        left_join(gene_biotypes, by = c('gene' = 'name')) |> 
        filter(biotype == 'lncRNA') |> 
        filter(chrom != 'mm39_X' & chrom != 'mm39_Y') |> 
        filter(abs(log2FoldChange) > 0.585) |> 
        filter(padj < 0.05)
    )
  } else if (list_type == "tyssowski") {
    PRG_rapid <- read_excel("02-data/published_data/Tyssowski2018/tyssowski_gene_lists.xlsx", sheet = 1) |> 
      pull(`Gene ID`)
    PRG_delay <- read_excel("02-data/published_data/Tyssowski2018/tyssowski_gene_lists.xlsx", sheet = 2) |> 
      pull(`Gene ID`)
    SRG <- read_excel("02-data/published_data/Tyssowski2018/tyssowski_gene_lists.xlsx", sheet = 3) |> 
      pull(`Gene ID`)
    
    gene_list <- c(PRG_rapid, PRG_delay)
    gene_list <- gene_list[!is.na(gene_list)] # remove NA
  } else if (list_type == 'DEGs') {
    # get all DEGs from the expression heatmap in figure 1
    gene_list <- read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/0_DEG_classifications.csv") |> 
      pull(gene)
  } else if (list_type == 'circadian') {
    gene_list <- known_circadian_genes <- c(
      'Per1', 'Per2', 'Per3', 'Clock', 'Bmal1', 'Cry1', 'Cry2', 'Nr1d1', 'Nr1d2'
    )
  }
  
  return(gene_list)
}


LoadActivityColors <- function(dataset = 'Dec2024', palette='colorblind') {
  if (palette == 'colorblind') {
    palette_colors <- c("#1E88E5", "#D81B60", "#FFC107", "#098154FF", "#63a45e", "#a8c66c")
  } else if (palette == 'grayscale') {
    # palette_colors <- c('gray90', 'gray75', 'gray60', 'gray45', 'gray30', 'gray15')
    palette_colors <- c('gray70', 'gray60', 'gray50', 'gray40', 'gray30', 'gray25')
    # palette_colors <- c('gray65', 'gray60', 'gray55', 'gray50', 'gray45', 'gray40')
  }

  if (dataset == 'May2024') {
    custom_colors <- c('SE'    = palette_colors[1], 
                       'EE30m' = palette_colors[2],
                       'EE6h'  = palette_colors[3],
                       'KA1h'  = palette_colors[4])
  } else if (dataset == 'Dec2023') {
    custom_colors <- c('SE' = palette_colors[1], 
                       '1h' = palette_colors[2],
                       '6h' = palette_colors[3],
                       'KA' = palette_colors[4])
  } else if (dataset == 'Dec2024') {
    custom_colors <- c('SE'    = palette_colors[1],
                       'EE30m' = palette_colors[2],
                       'EE6h'  = palette_colors[3],
                       'KA30m' = palette_colors[4],
                       'KA1h'  = palette_colors[6],
                       'KA6h'  = palette_colors[5]
    )
  }

  return(custom_colors)
}


LoadZTColors <- function(palette = 5) {
  if (palette == 1) {
    zt_colors <- c(
      'ZT0' = '#FADF7F',
      'ZT4' = '#D9B26F',
      'ZT12' = '#A69658',
      'ZT16' = '#795C5F'
    )
  } else if (palette == 2) {
    zt_colors <- c(
      'ZT0' = '#D9B26F',
      'ZT4' = '#F08080',
      'ZT12' = '#5A0001',
      'ZT16' = '#313628'
    )
  } else if (palette == 3) {
    zt_colors <- c(
      'ZT0' = '#92374D',
      'ZT4' = '#DDCA7D',
      'ZT12' = '#727D71',
      'ZT16' = '#082D0F'
    )
  } else if (palette == 4) {
    zt_colors <- c(
      'ZT0' = '#C2DFE3',
      'ZT4' = '#9DB4C0',
      'ZT12' = '#5C6B73',
      'ZT16' = '#253237'
    )
  } else if (palette == 5) {
    zt_colors <- c(
      'ZT0' = '#FFC49B',
      'ZT4' = '#FFEFD3',
      'ZT12' = '#C4D6E5',
      'ZT16' = '#ADB6C4'
    )
  }

  return(zt_colors)
}


LoadSexColors <- function() {
  sex_colors <- c('M' = '#E66100', 
                  'F' = '#5D3A9B',
                  'male' = '#E66100',
                  'female' = '#5D3A9B')
  
  return(sex_colors)
}


LoadSubclassesToUse <- function(seurat_obj=NULL, ascertainment = 'custom', as_gigaclasses = F, cell_cutoff = 150) {
# returns a list of subclass_names to use. 
# subclasses with fewer than cell_cutoff are excluded.
  attach_package_once("Seurat")

  if (ascertainment == 'custom') {
    subclass_list <- c(
      '016 CA1-ProS Glut',
      '025 CA2-FC-IG Glut',
      '017 CA3 Glut',
      '023 SUB-ProS Glut',
      '031 CT SUB Glut',
      '033 NP SUB Glut',
      '037 DG Glut',
      
      '053 Sst Gaba',
      '052 Pvalb Gaba',
      '051 Pvalb chandelier Gaba',
      '046 Vip Gaba',
      '047 Sncg Gaba',
      '049 Lamp5 Gaba',
      '050 Lamp5 Lhx6 Gaba',
      '048 RHP-COA Ndnf Gaba',
      
      '038 DG-PIR Ex IMN',
      '319 Astro-TE NN',
      '327 Oligo NN',
      '326 OPC NN',
      '334 Microglia NN'
    )
  } else if (ascertainment == 'auto') {
    subclass_list <- seurat_obj@meta.data |> 
      group_by(subclass_name) |> 
      summarize(n = n()) |> 
      filter(n >= cell_cutoff) |> 
      arrange(subclass_name) |> 
      pull(subclass_name)
  }
  
  if (as_gigaclasses) {
    subclass_sets <- c(
      excitatory = list(subclass_list[1:7]),
      inhibitory = list(subclass_list[8:15]),
      glia = list(subclass_list[c(16:20)])
    )
    
    return(subclass_sets)
  } else {
    return(subclass_list)
  }
}


LoadBarebonesTheme <- function(legend_position = 'none', ticks = 'none') {
  # Returns a minimal ggplot2 theme with as few elements as possible. Helpful for
  # Importing vector images into Illustrator.
  attach_package_once("ggplot2")
  
  # set tick parameters
  if (ticks == 'none') {
    ticks_param <- theme(axis.ticks = element_blank(), 
                         axis.line = element_blank())
  } else if (ticks == 'both') {
    ticks_param <- theme(axis.ticks = element_line(), 
                         axis.line = element_line())
  } else if (ticks == 'x') {
    ticks_param <- theme(axis.ticks.x = element_line(),
                         axis.line.x = element_line(),
                         axis.ticks.y  = element_blank(),
                         axis.line.y = element_blank())
  } else if (ticks == 'y') {
    ticks_param <- theme(axis.ticks.y = element_line(), 
                         axis.line.y = element_line(),
                         axis.ticks.x  = element_blank(),
                         axis.line.x = element_blank())
  }
  
  theme_barebones <- theme_classic() +
    theme(
      legend.position = legend_position,
      title = element_blank(),
      text = element_blank(),
      line = element_blank(),
      rect = element_blank(),
      panel.grid = element_blank(),
    ) +
    ticks_param
  
  return(theme_barebones)
}


GetSubclassContrasts <- function(seurat_obj, subclass, cell_cutoff = 30) {
# returns a list of activity_condition contrasts to use
# given they meet the cell_cutoff criterion
  attach_package_once("Seurat")

  acceptable_contrasts  <- c(
    'EE30m_SE',
    'EE6h_SE',
    'KA30m_SE',
    'KA6h_SE',
    'KA30m_EE30m',
    'KA6h_EE6h'
  )

  valid_conditions <- seurat_obj@meta.data |> 
    filter(subclass_name == subclass) |> 
    group_by(activity_condition) |> 
    summarize(n = n()) |> 
    filter(n >= 30) |> 
    pull(activity_condition)
  
  # make and filter contrast df
  df_contrasts <- expand.grid(group1 = valid_conditions, group2 = valid_conditions) |> 
    filter(group1 != group2) |> 
    mutate(pasted = glue('{group1}_{group2}')) |> 
    filter(pasted %in% acceptable_contrasts)
  
  return(df_contrasts)
}


GetCorrData <- function(seurat_obj, specific_condition, gene_list, output_fmt = 'complex_heatmap') {
# returns the pairwise correlation between expression values of all gene pairs
# in gene_list. Only cells in the specified condition will be included.
  attach_package_once("Seurat")

  expression_data <- seurat_obj |> 
    subset(activity_condition == specific_condition) |> 
    GetAssayData('SCT', layer = 'data') |> 
    as.data.frame() |>
    rownames_to_column(var = 'gene') |> 
    filter(gene %in% gene_list) |> 
    column_to_rownames(var = 'gene') |>
    as.matrix()

  if (output_fmt == 'complex_heatmap') {
    # format output as matrix
    corr <- cor(t(expression_data), method = 'pearson', use = 'pairwise.complete.obs')
    
    # find rows and columns that are exclusively NA and remove them
    corr <- corr[rowSums(is.na(corr)) < ncol(corr), ]
    corr <- corr[, colSums(is.na(corr)) < nrow(corr)]
    
    
    
  } else if (output_fmt == 'geom_tile') {
    # format output as df
      corr <- cor(t(expression_data), method = 'pearson') |>
        reshape2::melt() |>
        rename('gene1' = 'Var1', 'gene2' = 'Var2', 'corr_coef' = 'value') |>
        as_tibble() |>
        mutate(gene1 = factor(gene1, levels = gene_list),
               gene2 = factor(gene2, levels = gene_list)) |>
        arrange(gene1, gene2) |>
        mutate(condition = specific_condition) |>
        print()
  }
  
  return(corr)
}


PlotComplexHeatmap <- function(corr_matrix, plot_title, show_plot = TRUE, save_path = NULL) {
# Plots complex heatmap of genes x gene co-expression correlation values.
# Depends on output from GetCorrData, which must be called with argument
# output_fmt = 'complex_heatmap'.

  attach_package_once("Seurat")
  attach_package_once("ComplexHeatmap")
  attach_package_once("circlize")

  if (!is.null(save_path)) {
    # if (str_sub(save_path, -1) != '/') {
    #   save_path <- paste0(save_path, '/')
    # }
    # fname <- str_replace_all(plot_title, ' ', '_')
    # full_path <- paste0(save_path, fname, '.png')
    # 
    png(file = save_path, width = 1000, height = 900, units = 'px')
  }
  
  # make IEG annotation
  ieg_symbols <- LoadGeneList('IEG')
  bottom_anno <- HeatmapAnnotation(
    IEG = anno_mark(
      at = which(rownames(corr_matrix) %in% ieg_symbols),
      labels = intersect(rownames(corr_matrix), ieg_symbols),
      side = 'bottom'      
    ),
    show_annotation_name = FALSE
  )
  right_anno <- HeatmapAnnotation(
    IEG = anno_mark(
      at = which(colnames(corr_matrix) %in% ieg_symbols),
      labels = intersect(colnames(corr_matrix), ieg_symbols),
      side = 'bottom'      
    ),
    show_annotation_name = FALSE,
    which = 'row'
  )
  
  # plot using ComplexHeatmap library
  hm <- Heatmap(
    corr_matrix,
    bottom_annotation = bottom_anno,
    right_annotation = right_anno,
    column_title = plot_title,
    row_names_gp = gpar(fontsize = 20),
    column_names_gp = gpar(fontsize = 20),
    show_row_names = FALSE,
    show_column_names = FALSE,
    col = colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))
  )
  
  draw(hm)
  
  if (!is.null(save_path)) {dev.off()}

  return(hm)
}


FindActiveCells <- function(seurat_obj, subclass = '016 CA1-ProS Glut', gene_list = LoadGneList('IEG'), gene_threshold = 3) {
# Returns a tibble of cells that are "active" according to the following
# criterion: >= gene_threshold IEGs in a cell are expressed at counts >= 90th
# percentile of expression in the standard-environment condition.
  # find the 90th percentile of IEG expression in the dSE cells
  activation_thresholds <- seurat_obj |> 
    subset(subclass_name == subclass) |> 
    subset(activity_condition == 'SE') |> 
    GetAssayData('SCT', layer = 'data') |> 
    as.data.frame() |> 
    rownames_to_column(var = 'gene') |> 
    filter(gene %in% gene_list) |> 
    rowwise() |> 
    mutate(percentile_90th = quantile(c_across(-gene), 0.9)) |> 
    dplyr::select(gene, percentile_90th)
  
  activation_thresholds <- setNames(activation_thresholds$percentile_90th, activation_thresholds$gene)
  
  cell_conditions <- seurat_obj |> 
    FetchData("activity_condition") |> 
    rownames_to_column(var = "cell")
  
  # make a new df with only the expression for genes in gene_list
  df_gene_expression <- seurat_obj |> 
    subset(subclass_name == subclass) |>
    GetAssayData('SCT', layer = 'data') |> 
    as.data.frame() |> 
    rownames_to_column(var = 'gene') |>
    filter(gene %in% gene_list) |> 
    pivot_longer(cols = -gene, names_to = 'cell', values_to = 'expression') |> 
    pivot_wider(names_from = gene, values_from = expression)
  
  # confirm whether colnames in df_gene_expression are the same order as 
  # elements in activation_thresholds. this is important for comparing 
  # whether IEG expression exceeds 90th percentile ascertained thru SE
  df_colnames <- df_gene_expression |> 
    dplyr::select(-cell) |>
    colnames()
  
  verification <- all(df_colnames == names(activation_thresholds))
  if (verification != TRUE) {
    stop('colnames in df are not the same order as elements in thresholds')
  }
  
  # find out how many cells have at least 3 genes above activation_thresholds$percentile_90th
  df_active_cells <- df_gene_expression |> 
    rowwise() |> 
    mutate(num_upregd_genes = sum(c_across(-cell) > activation_thresholds)) |>
    left_join(cell_conditions, by = 'cell') |>
    mutate(active_binary = num_upregd_genes >= gene_threshold) |>
    relocate(activity_condition, active_binary, num_upregd_genes, .after = cell)
  
  # print percentages for each condition
  df_percent_active <- df_active_cells |>
    group_by(activity_condition) |> 
    summarize(percent_active = sum(active_binary) / n()) |> 
    ungroup() |> 
    mutate(subclass = subclass) |> 
    relocate(subclass)
 
  # package outputs
  results <- list(
    df_active_cells = df_active_cells,
    df_percent_active = df_percent_active,
    activation_thresholds = activation_thresholds
  )
  
  return(results)
}


FindDEGs <- function(seurat_obj, subclass, ident_var, group1, group2, logFC_threshold=0.25) {
# uses MAST to find DEGs between conditions within a subclass and removes sex-specific DEGs
  attach_package_once("Seurat")
  attach_package_once("glue")
  
  # print report statement
  print(glue("Finding DEGs between {subclass} {group1} and {group2}."))
  
  # subset nuclei
  nuclei_subclass <- seurat_obj |> 
    subset(subclass_name == subclass)
  
  # change cell identity to seurat_obj$ident_var
  ident_var_vector <- nuclei_subclass@meta.data |> pull(ident_var)
  Idents(nuclei_subclass) <- ident_var_vector
  
  # find DEGs
  DEGs <- FindMarkers(
    nuclei_subclass, 
    logfc.threshold = logFC_threshold,
    test.use = 'MAST',
    ident.1 = group1, ident.2 = group2,
    latent.vars = c('sublibrary', 'percent.mt'),
  ) |> 
    filter(p_val_adj < 0.05) |> 
    arrange(desc(avg_log2FC)) |> 
    rownames_to_column(var = 'gene')
  
  return(DEGs)
}


RunGOEnrichment <- function(target_gene_symbols, background_gene_list) {
  # Finds enriched GO terms of a target gene list run against the background genes.
  # target_gene_symbols should be, as implied, gene symbols. They are converted to 
  # Ensembl IDs enrichGO is run. Likewise, background_gene_list should also be in
  # symbol annotation.
  attach_package_once("clusterProfiler")
  attach_package_once("org.Mm.eg.db")
  attach_package_once("dplyr")

  # verify input is a gene list
  if (!is.character(target_gene_symbols)) {
    stop('input target_gene_symbols for RunGOEnrichment() must be a character vector.')
  } else if (!is.character(background_gene_list)) {
    stop('input background_gene_list for RunGOEnrichment() must be a character vector.')
  }
  
  # convert gene symbols to Ensembl IDs
  target_gene_IDs     <- SymbolToEnsembl(gene_symbol_vector = target_gene_symbols)
  background_gene_IDs <- SymbolToEnsembl(gene_symbol_vector = background_gene_list)
    
  # run analysis
  results_GO <- enrichGO(gene         = target_gene_IDs,
                  universe     = background_gene_IDs,
                  OrgDb        = org.Mm.eg.db,
                  keyType      = "ENSEMBL",
                  ont          = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.2) |> 
    setReadable(OrgDb = org.Mm.eg.db) |> 
    arrange(p.adjust) |> 
    print()
  
  return(results_GO@result)
}


RunYaoDGE <- function(input_subclass, save_to_file = FALSE) {
  attach_package_once("Seurat")
  
  # if yao not in namespace, load object
  if (!exists('yao')) {
    yao <- LoadDataset('Yao2023')
  }
  
  subclass_str <- gsub(' ', '_', input_subclass)
  subclass_str <- gsub('/', '', subclass_str)
  fname <- paste0("/Volumes/jack/yao/DGE/", subclass_str, "_DEGs.csv")
  
  if (file.exists(fname)) {
    return(read_csv(fname))
  }
  
  # subset subclass
  subclass_cells <- yao |> subset(subclass == input_subclass)
  
  # remove duplicate rownames
  gs <- rownames(ca1)
  idx <- rownames(ca1) |> duplicated() |> which()
  duplicate_gene_symbols <- gs[idx]
  
  # remove genes whose symbols are duplicated
  raw_counts <- GetAssayData(subclass_cells, assay = 'RNA')
  raw_counts <- raw_counts[!rownames(raw_counts) %in% duplicate_gene_symbols,]
  raw_counts <- CreateAssayObject(raw_counts, assay = 'RNA')
  
  # re-compose the seurat object
  subclass_meta <- subclass_cells@meta.data
  subclass_cells <- CreateSeuratObject(counts = raw_counts, meta.data = subclass_meta)
  
  # run DGE with MAST
  Idents(subclass_cells) <- subclass_cells$donor_sex
  
  DEGs <- subclass_cells |> 
    FindMarkers(ident.1 = 'F', ident.2 = 'M', 
                min.pct = 0.05,
                logfc.threshold = 0.1, 
                test.use = 'MAST',
                latent.vars = 'library_label') |>
    rownames_to_column(var = 'gene') |>
    filter(p_val_adj < 0.05)
  
  # save results as csv
  if (save_to_file) {
    write_csv(DEGs, fname)
  }
  
  return(DEGs)
}


QuickPercentExpression <- function(seurat.object, genes, ident_var, ident.1, ident.2, assay = "SCT", layer = "data") {
  # set Idents
  Idents(seurat.object) <- seurat.object[[ident_var]] |> pull()
  
  # Extract the available genes in the specified assay and layer
  available_genes <- rownames(GetAssayData(seurat.object, assay = assay, layer = layer))
  
  # Check if any of the requested genes are missing
  missing_genes <- setdiff(genes, available_genes)
  if(length(missing_genes) > 0) {
    warning("The following genes were not found in the assay data and will be skipped: ", 
            paste(missing_genes, collapse = ", "))
    genes <- intersect(genes, available_genes)
  }
  
  # Extract the expression matrix for the selected genes
  expr <- GetAssayData(seurat.object, assay = assay, layer = layer)[genes, , drop = FALSE]
  
  # Identify the cells for each group based on the provided identities
  cells1 <- WhichCells(seurat.object, idents = ident.1)
  cells2 <- WhichCells(seurat.object, idents = ident.2)
  
  # Compute the percentage of cells with expression > 0 in each group
  pct.1 <- rowSums(expr[, cells1] > 0) / length(cells1)
  pct.2 <- rowSums(expr[, cells2] > 0) / length(cells2)
  pct.mean <- rowMeans(cbind(pct.1, pct.2))
  shrink_coef <- sqrt(sqrt(pct.mean))
  
  # Create and return a data frame with the results
  result <- tibble(gene = genes, 
                   subclass = subclass,
                   group1 = ident.1,
                   group2 = ident.2,
                   pct.1 = pct.1, 
                   pct.2 = pct.2, 
                   pct.mean = pct.mean,
                   shrink_coef = shrink_coef)
  return(result)
}


PrintScriptDone <- function() {
  attach_package_once("glue")

  done_string <- "print(glue('Script {basename(sys.frame(1)$ofile)} complete!'))"
  return(done_string)  # run with eval(parse(done_string))
}


ShrinkSubclassName <- function(subclass_name) {
# removes slashes and underscores from Allen taxonomy subclass names 
# to prep for storage in filenames
  subclass_name_new <- gsub(" ", "_", subclass_name)
  subclass_name_new <- gsub("/", "", subclass_name_new)
  
  return(subclass_name_new)
}



MakeVolcanoPlot <- function(df, plot_title='log2(FC) vs -log10(p_val_adj)') {
  attach_package_once("dplyr")
  attach_package_once("plotly")
  attach_package_once("ggplot2")
  attach_package_once("ggrepel")
  
  # make static volcano plot first with ggplot
  threshold_line <- -log10(0.05)
  p <- ggplot(df, aes(x = avg_log2FC, y = -log10(p_val_adj), color = p_val_adj < 0.05)) +
    geom_point() +
    geom_text_repel(aes(label = gene)) +
    geom_hline(yintercept = threshold_line, linetype = 'dashed') +
    geom_vline(xintercept = 0) +
    # expand scale a little bit
    scale_x_continuous(expand = c(0.05, 0.05)) +
    scale_y_continuous(expand = c(0.05, 0.05)) +
    scale_color_manual(values = c('black', 'red')) +
    labs(title = plot_title, x = 'log2FC', y = '-log10(p_val_adj)') +
    theme_minimal()
  print(p)
  
  # Define colors explicitly to avoid RColorBrewer warnings
  df <- df |> mutate(color = ifelse(p_val_adj < 0.05, 'red', 'black'))
  
  # Create a formatted text column including multiple variables
  df <- df |> mutate(tooltip_text = paste0(
    "Gene: ", gene, "<br>",
    "log2FC: ", round(avg_log2FC, 3), "<br>",
    "pct.1: ", round(pct.1, 3), "<br>",
    "pct.2: ", round(pct.2, 3), "<br>",
    "adj. p-value: ", p_val_adj
  ))
  
  plot_ly(df, x = ~avg_log2FC, y = ~-log10(p_val_adj), 
          text = ~tooltip_text, hoverinfo = "text", 
          color = ~color, colors = c("black", "red")) |> 
    add_markers() |>
    layout(
      title = plot_title,
      xaxis = list(title = 'log2FC'),
      yaxis = list(title = '-log10(p_val_adj)'),
      shapes = list(
        list(
          type = "line",
          x0 = min(df$avg_log2FC, na.rm = TRUE), 
          x1 = max(df$avg_log2FC, na.rm = TRUE), 
          y0 = threshold_line, 
          y1 = threshold_line,
          line = list(color = "black", width = 1, dash = "dash")
        )
      )
    )
}



# GrowSubclassName <- function(subclass_name) {
#   # removes slashes and underscores from Allen taxonomy subclass names to prep for storage in filenames
#   subclass_name_new <- gsub(" ", "_", subclass_name)
#   subclass_name_new <- gsub("/", "", subclass_name_new)
#   
#   return(subclass_name_new)
# }


#' Patch of SeuratDisk::SaveH5Seurat function
#'
#' The "Assay5" class attribute "RNA" needs to be converted to a standard "Assay"
#' class for compatibility with SeuratDisk. It requires to make a temporary copy
#' so the size of the object grows bigger.
#'
#' @param object the Seurat object
#' @param filename the file path where to save the Seurat object
#' @param verbose SaveH5Seurat verbosity
#' @param overwrite whether to overwrite an existing file
#' 
#' @return NULL
SaveH5SeuratObject <- function(
    object,
    filename,
    verbose = TRUE,
    overwrite = TRUE
) {
  attach_package_once("Seurat")
  attach_package_once("SeuratDisk")

  # add copy of "RNA" 
  object[["RNA.tmp"]] <- CreateAssayObject(counts = object[["RNA"]]$counts)
  # switch default assay
  DefaultAssay(object) <- "RNA.tmp"
  # remove original
  object[["RNA"]] <- NULL
  # export
  SaveH5Seurat(object, filename = filename, overwrite, verbose)
  
  return(NULL)
}


#' Patch of SeuratDisk::LoadH5Seurat function
#'
#' The "Assay" class attribute "RNA" needs to be converted to the new "Assay5"
#' class. It requires to make a temporary copy so the size of the object grows
#' bigger.
#'
#' @param filename the file path where to save the Seurat object
#' @param verbose LoadH5Seurat verbosity
#' 
#' @return NULL
#' 
LoadH5SeuratObject <- function(filename, verbose = TRUE) {
  attach_package_once("Seurat")
  attach_package_once("SeuratDisk")

  # load h5 data
  object <- LoadH5Seurat(filename)
  # create "Assay5" class from old "Assay" class
  keys <- Key(object)
  slotID <- names(keys)[startsWith(keys, "rna")]
  object[["RNA"]] <- CreateAssay5Object(counts = object[[slotID]]$counts)
  # switch default assay
  DefaultAssay(object) <- "RNA"
  # delete "Assay" class
  object[[slotID]] <- NULL
  # reorder assays list
  object@assays <- object@assays[sort(names(object@assays))]
  
  return(object)
}


svgsave <- function(
    plot, 
    filename = 'tmp.svg', 
    path = '~/Downloads/', 
    height = 5, width = 5, 
    bg_color = 'transparent'
) {
  # Saves the current plot as an SVG file with the specified dimensions.
  attach_package_once("svglite")
  
  # clean filename by checking for .svg extension with regexp
  if (str_detect(filename, '.svg$') == F) {
    filename <- paste0(filename, '.svg')
  }
  # clean path
  if (str_sub(path, -1) != '/') {
    path <- paste0(path, '/')
  }
  full_path <- paste0(path, filename)
  
  svglite(filename = full_path, 
          width = width,
          height = height,
          bg = bg_color)
  print(plot)
  dev.off()
}


si <- function(width = 800, height = 800, format = 'svg', bg = 'white') {
  # Generate filename with an incrementing number
  fname <- "~/Downloads/tmp_1.svg"
  i <- 1
  while (file.exists(glue("~/Downloads/tmp_{i}.svg")) | file.exists(glue("~/Downloads/tmp_{i}.png"))) {
    i <- i + 1
    fname <- glue("~/Downloads/tmp_{i}.svg")
  }
  
  # toggle SVG vs PNG
  if (format == 'png') {
    fname <- gsub('svg', 'png', fname)
    dev.copy(png, file = fname, width = width, height = height, bg = bg)
    dev.off()
  } else {
    # Open SVG device, copy current plot, and close
    dev.copy(svg, file = fname, width = width / 100, height = height / 100, bg = bg)  # Convert pixels to inches
    dev.off()
  }
}
