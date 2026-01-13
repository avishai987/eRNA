#!/bin/bash
REGULATORY_GFF=$1        # GFF3 file with regulatory elements (e.g. enhancers)
OUTPUT_BED=$2      # Output BED file

# Convert GFF3 regulatory elements to BED format, keeping only enhancers

awk -v OFS="\t" '$3 == "enhancer" {
    # Try to extract the ID from the 9th column (e.g., ID=ENSR1_53X2;color=...)
    id = "NA"
    if (match($9, /ID=([^;]+)/, arr)) {
        id = arr[1]
    }
    print $1, $4-1, $5, id, $3
}' "$REGULATORY_GFF" > "$OUTPUT_BED"