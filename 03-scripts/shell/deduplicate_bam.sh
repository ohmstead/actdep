input_bam='/Volumes/jack/seq/analysis_May2024/sublibrary1/process/barcode_headAligned_anno_sorted_sublib1_cutoff80.bam'
output_bam=${input_bam%.bam}_dedup.bam

input_dir=$(dirname $input_bam)
log_file=$output_bam.log

echo "Input BAM file: $(basename $input_bam)"

# deduplicate the BAM file
echo "Deduplicating the BAM file..."
umi_tools dedup -I $input_bam -S $output_bam \
    --extract-umi-method=tag \
    --umi-tag=pN \
    --cell-tag=CB \
    --mapping-quality=30 2>&1 | tee $log_file

# index the output BAM file
echo "Indexing the output BAM file..."
samtools index $output_bam

echo "Done!"