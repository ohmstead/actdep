#!/bin/bash

input_bam="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/barcode_headAligned_anno_sorted_sublib1_cutoff80_dedup.bam"
barcode_list="04-analysis/May2024_barcodes.txt"
output_dir="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/sc_bams"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Count the total number of barcodes
total_barcodes=$(wc -l < "$barcode_list")
counter=0

# Loop through each barcode in the barcode list
while IFS= read -r barcode; do
    printf "Processing barcode %s (%d/%d)\n" "$barcode" "$((counter + 1))" "$total_barcodes"
    samtools view -h -d CB:$barcode "$input_bam" -@ 10 > "$output_dir/$barcode.bam"
    counter=$((counter + 1))
done < "$barcode_list"