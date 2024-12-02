#!/bin/bash

split-pipe \
--mode all \
--chemistry v3 \
--genome_dir /newvol/genomes/mm39 \
--fq1 /newvol/rawdata/JO1_S1_L006_R1_001.fastq.gz \
--output_dir /newvol/analysis/sublibrary1 \
--sample dSE1 A1 \
--sample d30m1 A2 \
--sample d6h1 A3 \
--sample KA1 A4 \
--sample dSE2 A5 \
--sample d30m2 A6 \
--sample d6h2 A7 \
--sample KA2 A8 \
--sample dSE3 A9 \
--sample d30m3 A10 \
--sample d6h3 A11 \
--sample KA3 A12

split-pipe \
--mode all \
--chemistry v3 \
--genome_dir /newvol/genomes/mm39 \
--fq1 /newvol/rawdata/JO2_S2_L006_R1_001.fastq.gz \
--output_dir /newvol/analysis/sublibrary2 \
--sample dSE1 A1 \
--sample d30m1 A2 \
--sample d6h1 A3 \
--sample KA1 A4 \
--sample dSE2 A5 \
--sample d30m2 A6 \
--sample d6h2 A7 \
--sample KA2 A8 \
--sample dSE3 A9 \
--sample d30m3 A10 \
--sample d6h3 A11 \
--sample KA3 A12

split-pipe \
--mode comb \
--sublibraries /newvol/analysis/sublibrary1 \
               /newvol/analysis/sublibrary2 \
--output_dir /newvol/analysis/combined