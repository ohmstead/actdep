input_bam="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/baby_bam.bam"
temp_bam="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/temp_bam.bam"
output_dir="/Volumes/jack/seq/analysis_May2024/sublibrary1/process/sample_bams"

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# split the bam file into individual samples from bam files
echo "Splitting temporary bam file into bam files for each sample"
samtools split -d pS -f "%!.bam" -v -@8 $input_bam

echo "done!"