# I would like to write a script that compares the May2024 and Dec2024 ZT12 timepoints.
# In principle, these should be very similar datasets. I will compare the DEGs between each condition and subclass.

# load libs & data ----
library(patchwork)
library(tidyverse)
library(venndir)
library(glue)

# load DEGs ----
May2024_deg_dir <- "04-analysis/DEGs/May2024_activity_condition"
Dec2024_deg_dir <- "04-analysis/DEGs/Dec2024_pilot_ZT12_activity_condition"

# among both directories, find shared files
May2024_files <- list.files(May2024_deg_dir)
Dec2024_files <- list.files(Dec2024_deg_dir)
contrast_list <- c('EE30m_vs_SE', 'EE6h_vs_SE', 'KA30m_vs_SE', 'KA6h_vs_SE', 'KA1h_vs_SE')

# end goal: for each subclass, plot all the contrasts that they have in common in a patchwork
# load in May2024 ----
degs_May2024 <- tibble()
for (full_fname in May2024_files) {
  fname <- str_remove(full_fname, ".csv")
  # Extract the contrast str (last two underscore-separated parts in filename)
  contrast <- str_extract(fname, "[^_]+_[^_]+_[^_]+$")
  # Extract the subclass str (everything before the contrast in filename)
  subclass <- str_remove(fname, paste0("_", contrast, "$"))  
  
  # load data
  degs_tmp <- read_csv(glue("{May2024_deg_dir}/{full_fname}")) |> 
    mutate(avg_log2FC = as.numeric(avg_log2FC)) |>
    filter(abs(avg_log2FC) > 0.585) |> 
    mutate(subclass = subclass, contrast = contrast) |> 
    mutate(subclass_by_contrast = glue("{subclass} x {contrast}")) |> 
    mutate(contrast_by_dataset = glue("{contrast} x May2024")) |>
    mutate(dataset = 'May2024') |> 
    select(gene, subclass_by_contrast, dataset, subclass, contrast, avg_log2FC)
  
  degs_May2024 <- bind_rows(degs_May2024, degs_tmp)
}


# load in Dec2024 ----
degs_Dec2024 <- tibble()
for (full_fname in Dec2024_files) {
  fname <- str_remove(full_fname, ".csv")
  # Extract the contrast str (last two underscore-separated parts in filename)
  contrast <- str_extract(fname, "[^_]+_[^_]+_[^_]+$")
  # Extract the subclass str (everything before the contrast in filename)
  subclass <- str_remove(fname, paste0("_", contrast, "$"))  
  
  # load data
  degs_tmp <- read_csv(glue("{Dec2024_deg_dir}/{full_fname}")) |> 
    mutate(subclass = subclass, contrast = contrast) |> 
    mutate(subclass_by_contrast = glue("{subclass} x {contrast}")) |> 
    mutate(contrast_by_dataset = glue("{contrast} x Dec2024")) |>
    mutate(avg_log2FC = as.numeric(avg_log2FC)) |>
    mutate(dataset = 'Dec2024') |> 
    select(gene, subclass_by_contrast, dataset, subclass, contrast, avg_log2FC)
  
  degs_Dec2024 <- bind_rows(degs_Dec2024, degs_tmp)
}


# combine datasets ----
degs <- bind_rows(degs_May2024, degs_Dec2024) |> 
  filter(contrast %in% contrast_list)


# get shared subclasses ----
shared_subclasses <- degs |> 
  group_by(subclass) |>
  summarize(n_datasets = n_distinct(dataset)) |> 
  filter(n_datasets == 2) |> 
  pull(subclass)

# plot deg overlaps ----
# for each subclass, make the following venn diagrams
# 1. May2024 EE30m_vs_SE vs Dec2024 EE30m_vs_SE
# 2. May2024 EE6h_vs_SE  vs Dec2024 EE6h_vs_SE
# 3. May2024 KA1h_vs_SE  vs Dec2024 KA30m_vs_SE
# 4. May2024 KA1h_vs_SE  vs Dec2024 KA6h_vs_SE
contrasts_to_compare <- tibble(
  order = c(1,2,3,4),
  May2024 = c('EE30m_vs_SE', 'EE6h_vs_SE', 'KA1h_vs_SE',  'KA1h_vs_SE'),
  Dec2024 = c('EE30m_vs_SE', 'EE6h_vs_SE', 'KA30m_vs_SE', 'KA6h_vs_SE')
)

k <- 1  # subclass counter
n <- length(shared_subclasses)

# for each subclass, make the 4 venn diagrams
for (current_subclass in shared_subclasses) {
  # get the DEGs for this subclass
  degs_subclass <- degs |> 
    filter(subclass == current_subclass)
  
  # make the venn diagrams
  for (i in 1:nrow(contrasts_to_compare)) {
    print(glue("Making venn diagram {i}/4 for subclass {k}/{n}..."))
    
    # get the contrasts to compare
    current_contrast <- contrasts_to_compare[i,]
    
    # get the genes for each dataset
    genes_May2024 <- degs_subclass |> 
      filter(dataset == 'May2024') |> 
      filter(contrast == current_contrast$May2024) |>
      pull(gene)
    
    genes_Dec2024 <- degs_subclass |> 
      filter(dataset == 'Dec2024') |> 
      filter(contrast == current_contrast$Dec2024) |>
      pull(gene)
    
    # create strings for plot
    plot_title <- glue("**{current_subclass}**\n\n**{current_contrast$May2024}** vs **{current_contrast$Dec2024}**")
    plot_fname <- glue("{current_subclass}___{current_contrast$May2024}___{current_contrast$Dec2024}.png")
    plot_path  <- glue("04-analysis/compareDatasetDEGs/{plot_fname}")
    
    # make the venn diagram
    png(plot_path, width = 11, height = 8.5, units = 'in', res = 600)
    
    if (length(genes_May2024) == 0 || length(genes_Dec2024) == 0) {
      # Create a plot showing one set is empty
      png(plot_path, width = 11, height = 8.5, units = 'in', res = 600)
      plot(1, type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", bty = "n", main = plot_title)
      text(1, 1, "One DEG list is empty, no overlap", cex = 1.5)
      dev.off()
    } else {
      # Regular Venn diagram
      png(plot_path, width = 11, height = 8.5, units = 'in', res = 600)
      venndir(
        setlist = list(
          May2024 = genes_May2024, 
          Dec2024_ZT12 = genes_Dec2024
        ),
        main = plot_title,
        proportional = TRUE,
        label_style = 'lite box',
        show_labels = 'NCPi',
        show_items = 'item',
      )
      dev.off()
    }
  }
  
  k <- k + 1
}
