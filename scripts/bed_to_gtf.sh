#!/bin/bash
source ~/.bashrc

# Script to convert a BED file of enhancers to GTF format.
# Assumes the BED file has at least 4 columns:
# 1. chromosome (seqname)
# 2. start (0-based)
# 3. end
# 4. name (enhancer identifier)
# And that the strand is always "."

# Input BED file path from the first command-line argument
input_bed="$1"
# Output GTF file path from the second command-line argument
output_gtf="$2"

# Check if both input and output file paths are provided
if [ -z "$input_bed" ] || [ -z "$output_gtf" ]; then
  echo "Usage: $0 <input_BED_file> <output_GTF_file>"
  exit 1
fi

# Define standard GTF fields
source="USER_DEFINED_ENHANCERS" # Source of the annotation
feature="regulatory_region"      # Feature type in GTF
score="."
strand="."
frame="."

# Read each line from the input BED file
while IFS=$'\t' read -r seqname start end name others; do
  # Use the name from the BED file (column 4) as the gene_id and transcript_id
  gene_id="$name"
  transcript_id="$name.1"

  # Construct the attributes string for the GTF file
  attributes="gene_id \"$gene_id\"; transcript_id \"$transcript_id\"; gene_biotype \"enhancer\"; transcript_biotype \"enhancer\";"

  # Construct the GTF line
  # Note: GTF start coordinate is 1-based, so we add 1 to the BED start.
  gtf_line="$seqname\t$source\t$feature\t$((start + 1))\t$end\t$score\t$strand\t$frame\t$attributes"

  # Write the GTF line to the output file
  printf "%b\n" "$gtf_line" >> "$output_gtf"

done < "$input_bed"

# Print a confirmation message
echo "Conversion from $input_bed to $output_gtf completed."