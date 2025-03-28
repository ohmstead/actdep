# Creates a list of 3 Seurat objects, 1 for each gigaclass.
# Should only need to be run once in order to save the objects to disk.
#  1. excitatory
#  2. inhibitory
#  3. glia

source("03-scripts/R/seq_functions.R")

library(Seurat)

# report gigaclasses
print('Saving the following gigaclass Seurat objects')
print('  excitatory:')
for (subclass in gigaclasses$excitatory) {print(glue('    {subclass}'))}
print('  inhibitory:')
for (subclass in gigaclasses$inhibitory) {print(glue('    {subclass}'))}
print('  glia:')
for (subclass in gigaclasses$glia) {print(glue('    {subclass}'))}

# load all nuclei
nuclei <- LoadDataset("Dec2024")
gigaclasses <- LoadSubclassesToUse(as_gigaclasses = TRUE)

# subset gigaclass objects
nuclei_excitatory <- nuclei |> subset(subclass_name %in% gigaclasses$excitatory)
nuclei_inhibitory <- nuclei |> subset(subclass_name %in% gigaclasses$inhibitory)
nuclei_glia       <- nuclei |> subset(subclass_name %in% gigaclasses$glia)

# save locally
save_dir <- "04-analysis/Seurats/Dec2024/"
SaveSeuratRds(nuclei_excitatory, file = paste0(save_dir, "seurat_excitatory.Rds"))
SaveSeuratRds(nuclei_inhibitory, file = paste0(save_dir, "seurat_inhibitory.Rds"))
SaveSeuratRds(nuclei_glia, file = paste0(save_dir, "seurat_glia.Rds"))
