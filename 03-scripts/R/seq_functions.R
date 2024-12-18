# A script housing a collection of commonly called functions in the actdep project.

RunYaoDGE <- function(input_subclass, save_to_file = FALSE) {
  library(Seurat)

  # if yao not in namespace, load object
  if (!exists('yao')) {
    yao <- LoadSeuratRds("/Volumes/jack/yao/yao_seurat.rds")
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

  if (dataset %in% c('yao', 'Yao')) {
    # try to load from local file
    if (path.exists()){
      local_file_loc <- "04-analysis/Seurats/Yao2023/seurat.Rds"
    } else {
      data_path <- "/Volumes/jack/yao/"
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
    custom_colors <- c('dSE' = palette_colors[2], 
                       'd30m' = palette_colors[1],
                       'd6h' = palette_colors[3],
                       'KA' = palette_colors[4])
  } else if (dataset == 'Dec2023') {
    custom_colors <- c('SE' = palette_colors[2], 
                       '1h' = palette_colors[1],
                       '6h' = palette_colors[3],
                       'KA' = palette_colors[4])
  } else if (dataset == 'Dec2024') {
    custom_colors <- c('SE' = palette_colors[2],
                       'EE_30m' = palette_colors[1],
                       'EE_6h' = palette_colors[3],
                       'KA_30m' = palette_colors[4],
                       'KA_1h' = palette_colors[5],
                       'KA_6h' = palette_colors[6]
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


LoadSubclassesToUse <- function(seurat_obj) {
# returns a list of subclass_names to use. 
# subclasses with fewer than 100 cells are excluded.
  library(Seurat)
  library(tidyverse)

  subclass_list <- seurat_obj@meta.data |> 
    filter(n() >= 100, .by = 'subclass_name') |> 
    distinct(subclass_name)
  
  rownames(subclass_list) <- NULL
  # sort list alphabetically
  subclass_list <- subclass_list$subclass_name |> sort()
  
  return(subclass_list)
}


LoadGeneList <- function(list_type = "IEG") {
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
  }
  
  return(gene_list)
}


GetCorrData <- function(seurat_obj, specific_condition, gene_list, output_fmt = 'complex_heatmap') {
# returns the pairwise correlation between expression values of all gene pairs
# in gene_list. Only cells in the specified condition will be included.
  library(Seurat)
  library(tidyverse)

  expression_data <- seurat_obj |> 
    subset(condition == specific_condition) |> 
    GetAssayData('SCT', layer = 'data') |> 
    as.data.frame() |>
    rownames_to_column(var = 'gene') |> 
    filter(gene %in% gene_list) |> 
    column_to_rownames(var = 'gene') |>
    as.matrix()

  if (output_fmt == 'complex_heatmap') {
    # format output as matrix
    corr <- cor(t(expression_data), method = 'pearson')
    
    # remove from t the rows and cols with names in the character vector u
    na_names <- names(which(is.na(corr[1,])))
    corr <- corr[!(rownames(corr) %in% na_names), !(colnames(corr) %in% na_names)]
    
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


PlotComplexHeatmap <- function(corr_matrix, plot_title, save_plot = FALSE) {
# Plots complex heatmap of genes x gene co-expression correlation values.
# Depends on output from GetCorrData, which must be called with argument
# output_fmt = 'complex_heatmap'.

  library(Seurat)
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)

  colormap <- colorRamp2(c(-1, 0, 1), c('blue', 'white', 'red'))
  
  if (save_plot) {
    folder_path <- '05-results/figure1/raw_R_plots/'
    fname <- paste0(folder_path, 'complexHeatmap_', plot_title, '.png')
    
    png(file = fname, width = 16.25, height = 15, units = 'in', res = 600)
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
  
  if (save_plot) {dev.off()}
}


FindActiveCells <- function(seurat_obj, gene_list, gene_threshold = 3) {
# Returns a tibble of cells that are "active" according to the following
# criterion: >= gene_threshold IEGs in a cell are expressed at counts >= 90th
# percentile of expression in the standard-environment condition.
  # find the 90th percentile of IEG expression in the dSE cells
  activation_thresholds <- seurat_obj |> 
    subset(subclass_name == '016 CA1-ProS Glut') |> 
    subset(condition == 'SE') |> 
    GetAssayData('SCT', layer = 'counts') |> 
    as.data.frame() |> 
    rownames_to_column(var = 'gene') |> 
    filter(gene %in% gene_list) |> 
    rowwise() |> 
    mutate(percentile_90th = quantile(c_across(-gene), 0.9)) |> 
    select(gene, percentile_90th)
  
  activation_thresholds <- setNames(activation_thresholds$percentile_90th, activation_thresholds$gene)
  
  cell_conditions <- seurat_obj |> 
    FetchData("condition") |> 
    rownames_to_column(var = "cell")
  
  # make a new df with only the expression for genes in gene_list
  df_ieg_expression <- seurat_obj |> 
    subset(subclass_name == '016 CA1-ProS Glut') |>
    GetAssayData('SCT', layer = 'counts') |> 
    as.data.frame() |> 
    rownames_to_column(var = 'gene') |>
    filter(gene %in% gene_list) |> 
    pivot_longer(cols = -gene, names_to = 'cell', values_to = 'expression') |> 
    pivot_wider(names_from = gene, values_from = expression)
  
  # confirm whether colnames in df_ieg_expression are the same order as 
  # elements in activation_thresholds. this is important for comparing 
  # whether IEG expression exceeds 90th percentile ascertained thru SE
  df_colnames <- df_ieg_expression |> 
    select(-cell) |>
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
    relocate(condition, active_binary, num_upregd_genes, .after = cell)
  
  # print percentages for each condition
  df_active_cells |>
    group_by(condition) |> 
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
  
  # from gene symbols, get Entrez gene IDs
  gene_info <- getBM(
    attributes = c("mgi_symbol", "entrezgene_id"),
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
