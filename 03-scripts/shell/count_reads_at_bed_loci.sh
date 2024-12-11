#!/bin/bash
# This script counts the number of reads at each locus in a BED file for each
# BAM file in a directory. It will process all BAMs in a folder, whether they
# are single-barcode or sample BAMs, and can count reads at any loci set that is
# defined in the BED file. The output of the script is a set of text files, one
# for each BAM file, with the counts at each locus.

bam_dir="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/sc_bams"
bed_file="02-data/published_data/eRNAbase/eRNA_merge_mouse_nonred_converted.bed"
output_dir="04-analysis/eRNA-counts/"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Get the total number of BAM files
total_bams=$(ls "$bam_dir"/*.sorted.bam | wc -l)
count=0

# Process each BAM file and append counts to the output file
for bam_file in "$bam_dir"/*.sorted.bam; do
    barcode=$(basename "$bam_file" .sorted.bam)
    printf "Processing barcode %s (%d/%d)\n" "$barcode" "$count" "$total_bams"
    
    # Count reads at BED loci (using only first 3 columns)
    cut -f1,2,3 "$bed_file" | \
    bedtools multicov -bams "$bam_file" -bed - | \
    awk 'OFS="\t" {print $1"__"$2"__"$3, $4}' > "04-analysis/eRNA-counts/${barcode}_counts.txt"
    
    count=$((count + 1))
    if [ $count -ge 5 ]; then
        break
    fi
done