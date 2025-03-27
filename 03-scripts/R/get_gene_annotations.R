# convert gtf to csv
# download gtf file from ensembl:
# https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz

gtf_data <- import("~/Downloads/Mus_musculus.GRCm39.113.gtf")
genes <- gtf_data[gtf_data$type == "gene"]
gene_df <- mcols(genes)[, c("gene_id", "gene_name", "gene_source", "gene_biotype")]

# remove duplicates (just in case)
gene_df <- unique(as.data.frame(gene_df))

# set ENSMUG genes with no name to their ENSMUSG
gene_df$gene_name[is.na(gene_df$gene_name)] <- gene_df$gene_id[is.na(gene_df$gene_name)]

# save
write.csv(gene_df, 
          "04-analysis/gene_biotypes.csv", 
          row.names = FALSE)
