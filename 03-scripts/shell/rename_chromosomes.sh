#!/bin/bash

# Input BAM file
INPUT_BAM="/Users/jack/code/seq/manuscript_code/figure_hyper_alignment/repex/KA2.bam"
OUTPUT_BAM="/Users/jack/code/seq/manuscript_code/figure_hyper_alignment/repex/KA2_chromosomes_renamed.bam"
NUM_CPU_CORES=10

# Extract the header
samtools view -@ $NUM_CPU_CORES -H $INPUT_BAM > header.sam
samtools view -H $INPUT_BAM > header.sam

# Replace the chromosome naming convention in the header
sed -i '' -e 's/SN:mm39_1/SN:chr1/g' \
          -e 's/SN:mm39_2/SN:chr2/g' \
          -e 's/SN:mm39_3/SN:chr3/g' \
          -e 's/SN:mm39_4/SN:chr4/g' \
          -e 's/SN:mm39_5/SN:chr5/g' \
          -e 's/SN:mm39_6/SN:chr6/g' \
          -e 's/SN:mm39_7/SN:chr7/g' \
          -e 's/SN:mm39_8/SN:chr8/g' \
          -e 's/SN:mm39_9/SN:chr9/g' \
          -e 's/SN:mm39_10/SN:chr10/g' \
          -e 's/SN:mm39_11/SN:chr11/g' \
          -e 's/SN:mm39_12/SN:chr12/g' \
          -e 's/SN:mm39_13/SN:chr13/g' \
          -e 's/SN:mm39_14/SN:chr14/g' \
          -e 's/SN:mm39_15/SN:chr15/g' \
          -e 's/SN:mm39_16/SN:chr16/g' \
          -e 's/SN:mm39_17/SN:chr17/g' \
          -e 's/SN:mm39_18/SN:chr18/g' \
          -e 's/SN:mm39_19/SN:chr19/g' \
          -e 's/SN:mm39_X/SN:chrX/g' \
          -e 's/SN:mm39_Y/SN:chrY/g' \
          -e 's/SN:mm39_MT/SN:chrM/g' \
          -e 's/SN:mm39_JH584299.1/SN:chr5_JH584299v1_random/g' \
          -e 's/SN:mm39_GL456233.2/SN:chrX_GL456233v2_random/g' \
          -e 's/SN:mm39_JH584301.1/SN:chrY_JH584301v1_random/g' \
          -e 's/SN:mm39_GL456211.1/SN:chr1_GL456211v1_random/g' \
          -e 's/SN:mm39_GL456221.1/SN:chr1_GL456221v1_random/g' \
          -e 's/SN:mm39_JH584297.1/SN:chr5_JH584297v1_random/g' \
          -e 's/SN:mm39_JH584296.1/SN:chr5_JH584296v1_random/g' \
          -e 's/SN:mm39_GL456354.1/SN:chr5_GL456354v1_random/g' \
          -e 's/SN:mm39_JH584298.1/SN:chr5_JH584298v1_random/g' \
          -e 's/SN:mm39_JH584300.1/SN:chrY_JH584300v1_random/g' \
          -e 's/SN:mm39_GL456219.1/SN:chr7_GL456219v1_random/g' \
          -e 's/SN:mm39_GL456210.1/SN:chr1_GL456210v1_random/g' \
          -e 's/SN:mm39_JH584303.1/SN:chrY_JH584303v1_random/g' \
          -e 's/SN:mm39_JH584302.1/SN:chrY_JH584302v1_random/g' \
          -e 's/SN:mm39_GL456212.1/SN:chr1_GL456212v1_random/g' \
          -e 's/SN:mm39_JH584304.1/SN:chrUn_JH584304v1/g' \
          -e 's/SN:mm39_GL456379.1/SN:chrUn_GL456379v1/g' \
          -e 's/SN:mm39_GL456366.1/SN:chrUn_GL456366v1/g' \
          -e 's/SN:mm39_GL456367.1/SN:chrUn_GL456367v1/g' \
          -e 's/SN:mm39_GL456239.1/SN:chr1_GL456239v1_random/g' \
          -e 's/SN:mm39_GL456383.1/SN:chrUn_GL456383v1/g' \
          -e 's/SN:mm39_GL456385.1/SN:chrUn_GL456385v1/g' \
          -e 's/SN:mm39_GL456360.1/SN:chrUn_GL456360v1/g' \
          -e 's/SN:mm39_GL456378.1/SN:chrUn_GL456378v1/g' \
          -e 's/SN:mm39_MU069435.1/SN:chrUn_MU069435v1/g' \
          -e 's/SN:mm39_GL456389.1/SN:chrUn_GL456389v1/g' \
          -e 's/SN:mm39_GL456372.1/SN:chrUn_GL456372v1/g' \
          -e 's/SN:mm39_GL456370.1/SN:chrUn_GL456370v1/g' \
          -e 's/SN:mm39_GL456381.1/SN:chrUn_GL456381v1/g' \
          -e 's/SN:mm39_GL456387.1/SN:chrUn_GL456387v1/g' \
          -e 's/SN:mm39_GL456390.1/SN:chrUn_GL456390v1/g' \
          -e 's/SN:mm39_GL456394.1/SN:chrUn_GL456394v1/g' \
          -e 's/SN:mm39_GL456392.1/SN:chrUn_GL456392v1/g' \
          -e 's/SN:mm39_GL456382.1/SN:chrUn_GL456382v1/g' \
          -e 's/SN:mm39_GL456359.1/SN:chrUn_GL456359v1/g' \
          -e 's/SN:mm39_GL456396.1/SN:chrUn_GL456396v1/g' \
          -e 's/SN:mm39_GL456368.1/SN:chrUn_GL456368v1/g' \
          -e 's/SN:mm39_MU069434.1/SN:chr1_MU069434v1_random/g' \
          -e 's/SN:mm39_JH584295.1/SN:chr4_JH584295v1_random/g' header.sam
# Extract the alignments without the header
samtools view -@ $NUM_CPU_CORES $INPUT_BAM | \
samtools view $INPUT_BAM | \
sed -e 's/\tmm39_1\t/\tchr1\t/g' \
    -e 's/\tmm39_2\t/\tchr2\t/g' \
    -e 's/\tmm39_3\t/\tchr3\t/g' \
    -e 's/\tmm39_4\t/\tchr4\t/g' \
    -e 's/\tmm39_5\t/\tchr5\t/g' \
    -e 's/\tmm39_6\t/\tchr6\t/g' \
    -e 's/\tmm39_7\t/\tchr7\t/g' \
    -e 's/\tmm39_8\t/\tchr8\t/g' \
    -e 's/\tmm39_9\t/\tchr9\t/g' \
    -e 's/\tmm39_10\t/\tchr10\t/g' \
    -e 's/\tmm39_11\t/\tchr11\t/g' \
    -e 's/\tmm39_12\t/\tchr12\t/g' \
    -e 's/\tmm39_13\t/\tchr13\t/g' \
    -e 's/\tmm39_14\t/\tchr14\t/g' \
    -e 's/\tmm39_15\t/\tchr15\t/g' \
    -e 's/\tmm39_16\t/\tchr16\t/g' \
    -e 's/\tmm39_17\t/\tchr17\t/g' \
    -e 's/\tmm39_18\t/\tchr18\t/g' \
    -e 's/\tmm39_19\t/\tchr19\t/g' \
    -e 's/\tmm39_X\t/\tchrX\t/g' \
    -e 's/\tmm39_Y\t/\tchrY\t/g' \
    -e 's/\tmm39_MT\t/\tchrM\t/g' \
    -e 's/\tmm39_JH584299.1\t/\tchr5_JH584299v1_random\t/g' \
    -e 's/\tmm39_GL456233.2\t/\tchrX_GL456233v2_random\t/g' \
    -e 's/\tmm39_JH584301.1\t/\tchrY_JH584301v1_random\t/g' \
    -e 's/\tmm39_GL456211.1\t/\tchr1_GL456211v1_random\t/g' \
    -e 's/\tmm39_GL456221.1\t/\tchr1_GL456221v1_random\t/g' \
    -e 's/\tmm39_JH584297.1\t/\tchr5_JH584297v1_random\t/g' \
    -e 's/\tmm39_JH584296.1\t/\tchr5_JH584296v1_random\t/g' \
    -e 's/\tmm39_GL456354.1\t/\tchr5_GL456354v1_random\t/g' \
    -e 's/\tmm39_JH584298.1\t/\tchr5_JH584298v1_random\t/g' \
    -e 's/\tmm39_JH584300.1\t/\tchrY_JH584300v1_random\t/g' \
    -e 's/\tmm39_GL456219.1\t/\tchr7_GL456219v1_random\t/g' \
    -e 's/\tmm39_GL456210.1\t/\tchr1_GL456210v1_random\t/g' \
    -e 's/\tmm39_JH584303.1\t/\tchrY_JH584303v1_random\t/g' \
    -e 's/\tmm39_JH584302.1\t/\tchrY_JH584302v1_random\t/g' \
    -e 's/\tmm39_GL456212.1\t/\tchr1_GL456212v1_random\t/g' \
    -e 's/\tmm39_JH584304.1\t/\tchrUn_JH584304v1\t/g' \
    -e 's/\tmm39_GL456379.1\t/\tchrUn_GL456379v1\t/g' \
    -e 's/\tmm39_GL456366.1\t/\tchrUn_GL456366v1\t/g' \
    -e 's/\tmm39_GL456367.1\t/\tchrUn_GL456367v1\t/g' \
    -e 's/\tmm39_GL456239.1\t/\tchr1_GL456239v1_random\t/g' \
    -e 's/\tmm39_GL456383.1\t/\tchrUn_GL456383v1\t/g' \
    -e 's/\tmm39_GL456385.1\t/\tchrUn_GL456385v1\t/g' \
    -e 's/\tmm39_GL456360.1\t/\tchrUn_GL456360v1\t/g' \
    -e 's/\tmm39_GL456378.1\t/\tchrUn_GL456378v1\t/g' \
    -e 's/\tmm39_MU069435.1\t/\tchrUn_MU069435v1\t/g' \
    -e 's/\tmm39_GL456389.1\t/\tchrUn_GL456389v1\t/g' \
    -e 's/\tmm39_GL456372.1\t/\tchrUn_GL456372v1\t/g' \
    -e 's/\tmm39_GL456370.1\t/\tchrUn_GL456370v1\t/g' \
    -e 's/\tmm39_GL456381.1\t/\tchrUn_GL456381v1\t/g' \
    -e 's/\tmm39_GL456387.1\t/\tchrUn_GL456387v1\t/g' \
    -e 's/\tmm39_GL456390.1\t/\tchrUn_GL456390v1\t/g' \
    -e 's/\tmm39_GL456394.1\t/\tchrUn_GL456394v1\t/g' \
    -e 's/\tmm39_GL456392.1\t/\tchrUn_GL456392v1\t/g' \
    -e 's/\tmm39_GL456382.1\t/\tchrUn_GL456382v1\t/g' \
    -e 's/\tmm39_GL456359.1\t/\tchrUn_GL456359v1\t/g' \
    -e 's/\tmm39_GL456396.1\t/\tchrUn_GL456396v1\t/g' \
    -e 's/\tmm39_GL456368.1\t/\tchrUn_GL456368v1\t/g' \
    -e 's/\tmm39_MU069434.1\t/\tchr1_MU069434v1_random\t/g' \
    -e 's/\tmm39_JH584295.1\t/\tchr4_JH584295v1_random\t/g' > alignments.sam
# Convert the modified alignments back to BAM format
samtools view -@ $NUM_CPU_CORES -b -T header.sam alignments.sam > alignments.bam
samtools view -b -T header.sam alignments.sam > alignments.bam
# Reassemble the BAM file with the modified header
samtools reheader -@ $NUM_CPU_CORES header.sam alignments.bam > $OUTPUT_BAM
samtools reheader header.sam alignments.bam > $OUTPUT_BAM

# Clean up temporary files
rm header.sam alignments.sam alignments.bam

echo "Chromosome naming convention has been updated in $OUTPUT_BAM"