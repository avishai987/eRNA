#!/bin/bash
source ~/.bashrc

set -e  # Exit immediately if any command fails

# Input files
peaks=$1 # narrow peaks from macs2
gff=$2  # GFF3-like file with columns: chrom, id, type, start, end, or BED with: chrom,star,end,type
name_col_idx=$3
output_file=$4

# Step 1 – Convert GFF-like file to BED format
# BED format: chrom, start (0-based), end, type
ext="${gff##*.}"  # extract extension (after the last dot)

if [[ "$ext" == "gff" || "$ext" == "gff3" ]]; then
    # convert GFF to BED
    awk 'BEGIN{OFS="\t"} {print $1, $4-1, $5, $3}' "$gff" > annotations.bed
else
    # assume input is already BED, just copy
    cp "$gff" annotations.bed
fi

# Step 2 – Intersect peaks with annotations
# -wa: write original peak (from -a)
# -wb: write overlapping annotation (from -b)
bedtools intersect -a "$peaks" -b annotations.bed -wa -wb > intersect_result.bed

# Count total unique peaks
total_peaks=$(cut -f1-3 $peaks | sort -u | wc -l)

# Write header line with total peaks
echo "total of $total_peaks peaks, intersected with:" > "$output_file"

# Step 3 Count how many unique peaks fall in each annotation type
awk -v OFS="\t" -v col="$name_col_idx" '{print $1":"$2"-"$3, $(col)}' intersect_result.bed \
  | sort | uniq \
  | awk '{count[$2]++} END {for (type in count) print count[type], type}' \
  | sort -nr >> "$output_file"

# Step 4 Count number of peaks that do not overlap any annotation
bedtools intersect -a "$peaks" -b annotations.bed -v \
  | awk 'BEGIN{count=0} {count++} END {print count, "unannotated"}' >> "$output_file"

# Step 5 Clean up intermediate files
rm -f annotations.bed intersect_result.bed

echo "✅ Done. Output written to: $output_file"
