# Load required libraries
library(Seurat)
library(dplyr)
library(tidyr)
library(tibble)

subclass_list <- LoadSubclassesToUse(nuclei)
nuclei <- LoadDataset('Dec2024')
nuclei <- subset(nuclei, subclass_name %in% subclass_list & 
                         activity_condition %in% c('SE', 'EE30m', 'EE6h'))

# ------------------------------
# Step 1: Compute global gene-level statistics
# ------------------------------

# Extract normalized (log-transformed) expression data
nuclei_tmp <- subset(seurat_subsets$excitatory, activity_condition %in% c('SE', 'EE30m', 'EE6h'))
data_matrix <- GetAssayData(nuclei_tmp, slot = "data")

# Compute global mean and variance for each gene
gene_stats <- tibble(
  gene = rownames(data_matrix),
  mean_expr = rowMeans(data_matrix),
  var_expr = matrixStats::rowVars(as.matrix(data_matrix))
)

# Create expression bins (using deciles, but with a fine resolution to allow unique breaks)
quantiles <- unique(quantile(gene_stats$mean_expr, probs = seq(0, 1, by = 0.001), na.rm = TRUE))
gene_stats <- gene_stats %>%
  mutate(mean_bin = cut(mean_expr, breaks = quantiles, include.lowest = TRUE))

# ------------------------------
# Step 2: Define IEGs and subset summary stats
# ------------------------------

# Define your IEG list (adjust as needed)
IEGs <- LoadGeneList('IEG')
ieg_stats <- gene_stats %>% filter(gene %in% IEGs)

# ------------------------------
# Step 3: Vectorized Matching Approaches
# ------------------------------

# We'll perform a cross join (Cartesian join) between IEGs and all candidate genes.
# This is achieved by joining on an empty key.

# Nearest Neighbor Matching:
nn_matches <- ieg_stats %>%
  rename(IEG = gene, IEG_mean = mean_expr, IEG_variance = var_expr, IEG_bin = mean_bin) %>%
  inner_join(
    gene_stats %>% 
      rename(matched_gene = gene,
             candidate_mean = mean_expr,
             candidate_variance = var_expr,
             candidate_bin = mean_bin),
    by = character()
  ) %>%
  filter(IEG != matched_gene) %>%  # Exclude self-matches
  mutate(distance = sqrt((IEG_mean - candidate_mean)^2 + (IEG_variance - candidate_variance)^2)) %>%
  group_by(IEG) %>%
  arrange(distance, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  mutate(match_method = "nearest",
         similarity_rank = row_number()) %>%
  select(IEG, IEG_mean, IEG_variance, matched_gene, match_method, similarity_rank) |> 
  print()

# Binning Strategy:
bin_matches <- ieg_stats %>%
  rename(IEG = gene, IEG_mean = mean_expr, IEG_variance = var_expr, IEG_bin = mean_bin) %>%
  inner_join(
    gene_stats %>% 
      rename(matched_gene = gene,
             candidate_mean = mean_expr,
             candidate_variance = var_expr,
             candidate_bin = mean_bin),
    by = character()
  ) %>%
  filter(IEG != matched_gene, IEG_bin == candidate_bin) %>%  # Limit to same bin
  mutate(distance = sqrt((IEG_mean - candidate_mean)^2 + (IEG_variance - candidate_variance)^2)) %>%
  group_by(IEG) %>%
  arrange(distance, .by_group = TRUE) %>%
  slice_head(n = 5) %>%
  mutate(match_method = "binning",
         similarity_rank = row_number()) %>%
  select(IEG, IEG_mean, IEG_variance, matched_gene, match_method, similarity_rank) |> 
  print()

# Combine matching results from both methods
matching_df <- bind_rows(nn_matches, bin_matches) |> 
  print()


# ------------------------------
for (gene in IEGs) {
# gene <- IEGs[1]
p <- VlnPlot(nuclei_tmp,
        features = c(gene, nn_matches |> filter(IEG==gene) |> pull(matched_gene)), 
        group.by = 'activity_condition')
print(p)
}

# genes were selected from the above plots because they had roughly similar expression, 
# but did not vary between activity conditions that was apparent by visual inspection
matched_genes_excitatory <- c(
  'Gm19744',
  'Lepr',
  'Cdc42ep3',
  'A830036E02Rik',
  '4732463B04Rik',
  'Klhl4',
  'Pstpip2',
  'Clhc1',
  '4930448C13Rik',
  'Ednra',
  'Gm35835',
  'Cfap70',
  'Arhgef6',
  'Ctnna3',
  'Adgra1'
)
names(matched_genes_excitatory) <- IEGs
matched_genes_excitatory
