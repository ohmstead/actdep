#!/bin/bash

input_bam="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/barcode_headAligned_anno_sorted.bam"
barcode_list="/Users/jack/Downloads/barcodes.txt"
output_dir="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/sc_bams"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Function to process a single barcode
process_barcode() {
    barcode=$1
    echo "Processing barcode: $barcode"
    samtools view -h "$input_bam" | grep -E "@|CB:Z:$barcode" | samtools view -bS - > "$output_dir/$barcode.bam"
}

export -f process_barcode
export input_bam
export output_dir

# Use GNU parallel to process barcodes concurrently
parallel -j 8 process_barcode ::: $(cat $barcode_list)
