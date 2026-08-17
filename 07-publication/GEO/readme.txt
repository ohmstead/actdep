The processed data files are as follows:

	* count_matrix.mtx : This is the raw barcode-gene count matrix of "unfiltered" nuclei from Parse Biosciences's Trailmaker pipeline.

	* all_genes.csv : A csv of Ensemble IDs and gene symbols for the mouse GRCm39 v113 genome.

	* cell_metadata.csv : Basic metadata for nuclei barcodes including biological replicate of origin (i.e. "sample"), gene count, etc.

	* seurat_trailmaker_raw.rds : An R variable file downloaded from Parse Biosciences's Trailmaker pipeline containing a basic, "filtered" Seurat v5 object of ~110,000 nuclei.

	* seurat_final.rds : An R variable of processed single nuclei with extensive metadata including ABC taxonomy assignment used for all downstream analysis. May be used to replicate analyses from the paper.