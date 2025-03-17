library(glue)
library(Seurat)

source("03-scripts/R/seq_functions.R")

# load Seurat and relevant subclasses
nuclei <- LoadDataset("Dec2024")
subclass_list <- LoadSubclassesToUse(nuclei, ascertainment = 'custom')

# load entire gene list for all subclasses
df_expression <- read_csv("04-analysis/DEGs/Dec2024_activity_condition_pseudobulk/0_DEG_classifications.csv") |> 
  print()

# run loop to generate gene-gene correlation matricies
plots <- list()
i <- 1
for (subclass in subclass_list) {
  print(glue("Processing {i}/{length(subclass_list)}: {subclass}"))
  print('Subsetting Seurat object...')
  nuclei_subclass <- subset(nuclei, subclass_name == subclass)
  
  subclass_fname <- ShrinkSubclassName(subclass)
  plot_path <- glue("05-results/figure3/raw_R_plots/{subclass_fname}_EE30m.png")
  
  print('Getting gene-gene correlations...')
  tmp <- GetCorrData(nuclei_subclass, 'EE30m', df_expression |> distinct(gene) |> pull())
  
  print('Plotting and saving...')
  p <- PlotComplexHeatmap(tmp, plot_title = glue("{subclass} EE30m cells"), save_path = plot_path)
  
  plots[[subclass]] <- p
  i <- i + 1
}

htShiny(plots[[19]])
