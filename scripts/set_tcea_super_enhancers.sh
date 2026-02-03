#!/bin/bash

# Define command-line arguments
params_dir=$1
input_hg19ToHg38_over_chain=$2
output_super_enhancers_bed=$3

wget https://bioinformatics.mdanderson.org/Supplements/Super_Enhancer/5_Super_enhancer_annotation/Annotation_eRNA_in_Super_enhancer_all.bed.gz -P "$params_dir"
gunzip "$params_dir"Annotation_eRNA_in_Super_enhancer_all.bed.gz
#remove first row (header)
sed -i '1d' "$params_dir"Annotation_eRNA_in_Super_enhancer_all.bed
liftOver "$params_dir"Annotation_eRNA_in_Super_enhancer_all.bed "$input_hg19ToHg38_over_chain" "$output_super_enhancers_bed" "$params_dir"unMapped_superEnhancers.bed
# keep only standard chromosomes in the first column
awk '$1 ~ /^chr[0-9XYM]{1,2}$/' "$output_super_enhancers_bed" > "$output_super_enhancers_bed".tmp
mv "$output_super_enhancers_bed".tmp "$output_super_enhancers_bed"
#download chromosome info file
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz -P "$params_dir"
gunzip "$params_dir"chromInfo.txt.gz
# add 150 bp flanks to each side (annotation file only includes core enhancer of 1bp)
bedtools slop -i "$output_super_enhancers_bed" -g <(awk '{print $1"\t"$2}' "$params_dir"chromInfo.txt) -b 150 > "$output_super_enhancers_bed".tmp
mv "$output_super_enhancers_bed".tmp "$output_super_enhancers_bed"
# remove "chr" from first column
awk '{$1 = gensub(/^chr/, "", 1, $1); print}' OFS="\t" "$output_super_enhancers_bed" > "$output_super_enhancers_bed".tmp
mv "$output_super_enhancers_bed".tmp "$output_super_enhancers_bed"
# intersect with PBMC enhancers
bedtools intersect -a "$output_super_enhancers_bed" -b Analysis/10X_PBMC/cell_type_enhancers/combined_enhancers_loci.bed  -wa -f 0.7 > "$output_super_enhancers_bed".tmp
mv "$output_super_enhancers_bed".tmp "$output_super_enhancers_bed"