#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status

# Input files
FEATURES_GFF=$1             # GFF3 file containing various genomic features
REGULATORY_GFF=$2        # GFF3 file with regulatory elements (e.g. enhancers)
OUTPUT_BED=$3      # Final output file with filtered enhancers

# List of feature types to exclude (RNA-producing or gene-associated regions)
EXCLUDE_FEATURES="mRNA|transcript|exon"


# Temporary files
TMP_FEATURES_BED="features_exclude.bed"
TMP_ENHANCERS_BED="enhancers.bed"

# Step 1: Convert GFF3 features to BED format and filter only the features to exclude
# GFF format: columns 1 = chrom, 4 = start, 5 = end, 3 = feature type
awk -v OFS="\t" -v pat="$EXCLUDE_FEATURES" '$3 ~ pat {print $1, $4-1, $5}' "$FEATURES_GFF" > "$TMP_FEATURES_BED"

# Step 2: Convert GFF3 regulatory elements to BED format, keeping only enhancers

awk -v OFS="\t" '$3 == "enhancer" {
    # Try to extract the ID from the 9th column (e.g., ID=ENSR1_53X2;color=...)
    id = "NA"
    if (match($9, /ID=([^;]+)/, arr)) {
        id = arr[1]
    }
    print $1, $4-1, $5, id, $3 
}' "$REGULATORY_GFF" > "$TMP_ENHANCERS_BED"

# Step 3: Use bedtools to find enhancers that do NOT overlap any of the excluded features
bedtools intersect -a "$TMP_ENHANCERS_BED" -b "$TMP_FEATURES_BED" -v > "$OUTPUT_BED"

# Step 4: Clean up temporary files
rm "$TMP_FEATURES_BED" "$TMP_ENHANCERS_BED"


# Final message
echo "Filtered enhancers saved to: $OUTPUT_BED"
