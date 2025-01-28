# A script housing a collection of commonly called functions in the actdep project.

RunYaoDGE <- function(input_subclass, save_to_file = FALSE) {
  library(Seurat)

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


LoadDataset <- function(dataset, sublibrary = "combined") {
  library(Seurat)
  library(glue)

  if (dataset %in% c('yao', 'Yao', 'yao2023', 'Yao2023')) {
    # try to load from local file
    if (file.exists("04-analysis/Seurats/Yao2023/seurat.Rds")){
      local_file_loc <- "04-analysis/Seurats/Yao2023/seurat.Rds"
    } else {
      data_path <- "/Volumes/jack/Yao2023/"
      seurat_obj <- LoadSeuratRds(paste0(data_path, "yao_seurat.rds"))
    }
  } else {
    # try to load from local file
    local_file_loc <- glue("04-analysis/Seurats/{dataset}/seurat.Rds")
    if (file.exists(local_file_loc)) {
      seurat_obj <- LoadSeuratRds(local_file_loc)
    } else { # load from SSD
      ssd_file_loc <- glue("/Volumes/jack/seq/analysis_{dataset}/{sublibrary}/0_all-sample/DGE_filtered/seurat.Rds")
      seurat_obj <- LoadSeuratRds(ssd_file_loc)
    }
  }
  return(seurat_obj)
}


LoadAllenColors <- function(clade = 'subclass') {
  library(tidyverse)

  allen_colors <- read_csv("02-data/published_data/allen_taxonomy_colors.csv")
  
  switch(
    clade,
    'class' = return(allen_colors |> distinct(class, class_color) |> deframe()),
    'subclass' = return(allen_colors |> distinct(subclass, subclass_color) |> deframe()),
    'supertype' = return(allen_colors |> distinct(supertype, supertype_color) |> deframe()),
    'cluster' = return(allen_colors |> distinct(cluster, cluster_color) |> deframe())
  )
}


LoadSexColors <- function() {
  sex_colors <- c('M' = '#E66100', 
                  'F' = '#5D3A9B',
                  'male' = '#E66100',
                  'female' = '#5D3A9B')
  
  return(sex_colors)
}


LoadActivityColors <- function(dataset = 'Dec2024', palette='colorblind') {
  if (palette == 'colorblind') {
    palette_colors <- c("#D81B60", "#1E88E5", "#FFC107", "#098154FF", "#63a45e", "#a8c66c")
  } else {
    # palette_colors <- paletteer::paletteer_d('ggthemes::wsj_colors6')
  }

  if (dataset == 'May2024') {
    custom_colors <- c('SE' = palette_colors[2], 
                       'EE30m' = palette_colors[1],
                       'EE6h' = palette_colors[3],
                       'KA1h' = palette_colors[4])
  } else if (dataset == 'Dec2023') {
    custom_colors <- c('SE' = palette_colors[2], 
                       '1h' = palette_colors[1],
                       '6h' = palette_colors[3],
                       'KA' = palette_colors[4])
  } else if (dataset == 'Dec2024') {
    custom_colors <- c('SE' = palette_colors[2],
                       'EE30m' = palette_colors[1],
                       'EE6h' = palette_colors[3],
                       'KA30m' = palette_colors[4],
                       'KA1h' = palette_colors[5],
                       'KA6h' = palette_colors[6]
    )
  }

  return(custom_colors)
}


LoadZTColors <- function(palette = 1) {
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


LoadSubclassesToUse <- function(seurat_obj, cell_cutoff = 150) {
# returns a list of subclass_names to use. 
# subclasses with fewer than 100 cells are excluded.
  library(Seurat)
  library(tidyverse)

  subclass_list <- seurat_obj@meta.data |> 
    group_by(subclass_name) |> 
    summarize(n = n()) |> 
    filter(n >= cell_cutoff)
  
  # sort list alphabetically
  subclass_list <- subclass_list |> 
    arrange(subclass_name)
  
  return(subclass_list)
}


LoadGeneList <- function(list_type = "IEG") {
  library(readxl)
  library(dplyr)
  
  if (list_type == "IEG") {
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
    gene_list <- c(
      ""
    )
  } else if (list_type == "tyssowski") {
    PRG_rapid <- read_excel("02-data/published_data/Tyssowski2018/tyssowski_gene_lists.xlsx", sheet = 1) |> 
      pull(`Gene ID`)
    PRG_delay <- read_excel("02-data/published_data/Tyssowski2018/tyssowski_gene_lists.xlsx", sheet = 2) |> 
      pull(`Gene ID`)
    
    gene_list <- c(PRG_rapid, PRG_delay)
  }
  
  return(gene_list)
}


GetCorrData <- function(seurat_obj, specific_condition, gene_list, output_fmt = 'complex_heatmap') {
# returns the pairwise correlation between expression values of all gene pairs
# in gene_list. Only cells in the specified condition will be included.
  library(Seurat)
  library(tidyverse)

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


PlotComplexHeatmap <- function(corr_matrix, plot_title, save_path = NULL) {
# Plots complex heatmap of genes x gene co-expression correlation values.
# Depends on output from GetCorrData, which must be called with argument
# output_fmt = 'complex_heatmap'.

  library(Seurat)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)

  colormap <- colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))
  
  if (!is.null(save_path)) {
    if (str_sub(save_path, -1) != '/') {
      save_path <- paste0(save_path, '/')
    }
    fname <- str_replace_all(plot_title, ' ', '_')
    full_path <- paste0(save_path, 'complexHeatmap_', fname, '.png')
    
    png(file = full_path, width = 1000, height = 900, units = 'px')
  }
  
  # plot using ComplexHeatmap library
  hm <- Heatmap(
    corr_matrix, 
    column_title = plot_title,
    row_names_gp = gpar(fontsize = 30),
    column_names_gp = gpar(fontsize = 30),
    col = colormap
  )
  
  draw(hm)
  print(hm)
  
  if (!is.null(save_path)) {dev.off()}
}


FindActiveCells <- function(seurat_obj, gene_list, gene_threshold = 3) {
# Returns a tibble of cells that are "active" according to the following
# criterion: >= gene_threshold IEGs in a cell are expressed at counts >= 90th
# percentile of expression in the standard-environment condition.
  # find the 90th percentile of IEG expression in the dSE cells
  activation_thresholds <- seurat_obj |> 
    subset(subclass_name == '016 CA1-ProS Glut') |> 
    subset(activity_condition == 'SE') |> 
    GetAssayData('RNA', layer = 'data') |> 
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
  df_ieg_expression <- seurat_obj |> 
    subset(subclass_name == '016 CA1-ProS Glut') |>
    GetAssayData('RNA', layer = 'counts') |> 
    as.data.frame() |> 
    rownames_to_column(var = 'gene') |>
    filter(gene %in% gene_list) |> 
    pivot_longer(cols = -gene, names_to = 'cell', values_to = 'expression') |> 
    pivot_wider(names_from = gene, values_from = expression)
  
  # confirm whether colnames in df_ieg_expression are the same order as 
  # elements in activation_thresholds. this is important for comparing 
  # whether IEG expression exceeds 90th percentile ascertained thru SE
  df_colnames <- df_ieg_expression |> 
    dplyr::select(-cell) |>
    colnames()
  
  verification <- all(df_colnames == names(activation_thresholds))
  if (verification != TRUE) {
    stop('colnames in df are not the same order as elements in thresholds')
  }
  
  # find out how many cells have at least 3 genes above activation_thresholds$percentile_90th
  df_active_cells <- df_ieg_expression |> 
    rowwise() |> 
    mutate(num_upregd_genes = sum(c_across(-cell) > activation_thresholds)) |>
    left_join(cell_conditions, by = 'cell') |>
    mutate(active_binary = num_upregd_genes >= gene_threshold) |>
    relocate(activity_condition, active_binary, num_upregd_genes, .after = cell)
  
  # print percentages for each condition
  df_active_cells |>
    group_by(activity_condition) |> 
    summarize(percent_active = sum(active_binary) / n()) |> 
    ungroup() |> 
    print()
  
  # package outputs
  results <- list(
    df_active_cells = df_active_cells,
    activation_thresholds = activation_thresholds
  )
  
  return(results)
}


FindDEGs <- function(seurat_obj, subclass, ident_var, group1, group2, logFC_threshold=0.25) {
# uses MAST to find DEGs between conditions within a subclass and removes sex-specific DEGs
  library(Seurat)
  library(glue)
  
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


RunGOEnrichment <- function(gene_list){
# accepts gene_list as input, returns df of enriched GO terms
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(biomaRt)
  library(tidyverse)

  # verify input is a gene list
  if (!is.character(gene_list)) {
    stop('gene_list input for RunGOEnrichment() must be a character vector.')
  }
  
  # make ensemble biomart object
  ensembl <- useEnsembl(biomart = "ensembl", 
                     dataset = "mmusculus_gene_ensembl")
  
  
  # from gene symbols, get Entrez gene IDs
  gene_info <- getBM(
    attributes = c("mgi_symbol", "entrezgene_id", "ensembl_gene_id"),
    filters = "mgi_symbol",
    values = gene_list,
    mart = ensembl
  )
  
  # run analysis
  results_GO <- enrichGO(gene = gene_info$entrezgene_id, OrgDb = 'org.Mm.eg.db', ont = 'ALL')
  results_GO <- as_tibble(setReadable(results_GO, OrgDb = 'org.Mm.eg.db', keyType = 'ENTREZID'))
  
  return(results_GO)  
}


PrintScriptDone <- function() {
  library(glue)

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
  library(Seurat)
  library(SeuratDisk)

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
  library(Seurat)
  library(SeuratDisk)

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
