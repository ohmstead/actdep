#!/bin/bash

# /newvolume/rawdata_Feb2024/

# JO1_S1_L007_R1_001.fastq.gz
# JO1_S1_L007_R2_001.fastq.gz  

# JO2_S2_L007_R1_001.fastq.gz  
# JO2_S2_L007_R2_001.fastq.gz  

# /newvolume/rawdata_Dec2023/

# JO1_S1_L008_R1_001.fastq.gz  
# JO1_S1_L008_R2_001.fastq.gz

# JO2_S2_L008_R1_001.fastq.gz 
# JO2_S2_L008_R2_001.fastq.gz


# cat /newvolume/rawdata_Dec2023/JO1_S1_L007_R1_001.fastq.gz /newvolume/rawdata_Feb2024/JO1_S1_L008_R1_001.fastq.gz > /newvolume/rawdata_merged/JO1_R1_merge.fastq.gz
# cat /newvolume/rawdata_Dec2023/JO1_S1_L007_R2_001.fastq.gz /newvolume/rawdata_Feb2024/JO1_S1_L008_R2_001.fastq.gz > /newvolume/rawdata_merged/JO1_R2_merge.fastq.gz
# cat /newvolume/rawdata_Dec2023/JO2_S2_L007_R1_001.fastq.gz /newvolume/rawdata_Feb2024/JO2_S2_L008_R1_001.fastq.gz > /newvolume/rawdata_merged/JO2_R1_merge.fastq.gz
# cat /newvolume/rawdata_Dec2023/JO2_S2_L007_R2_001.fastq.gz /newvolume/rawdata_Feb2024/JO2_S2_L008_R2_001.fastq.gz > /newvolume/rawdata_merged/JO2_R2_merge.fastq.gz

cat /Volumes/jack/seq/rawdata_JO1/JO1_S1_L007_R1_001.fastq.gz /Volumes/jack/seq/rawdata_JO1/JO1_S1_L008_R1_001.fastq.gz > /Volumes/jack/seq/rawdata_JO1/JO1_R1_merge.fastq.gz
cat /Volumes/jack/seq/rawdata_JO1/JO1_S1_L007_R2_001.fastq.gz /Volumes/jack/seq/rawdata_JO1/JO1_S1_L008_R2_001.fastq.gz > /Volumes/jack/seq/rawdata_JO1/JO1_R2_merge.fastq.gz