#!/bin/bash
# Convert chromatin counts (all rows except the first) first column to bed format
# Usage: ./chromatin_regions_to_bed.sh <chromatin_file> <output_bed_file>
# remove all "chr" from chromosome names from column 1
#input example: chr2:148881654-148881927
#output example: 2	148881654	148881927	chr2:148881654-148881927
chromatin_file=$1
output_file=$2
cat "$chromatin_file" | tail -n +2 | cut -f1 | awk -F'[:-]' 'BEGIN {OFS="\t"} {chr=$1; gsub(/^chr/, "", chr); print chr, $2, $3, $0}' > "$output_file"