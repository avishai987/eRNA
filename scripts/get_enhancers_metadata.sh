#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <enhancers_bed> <gff3_file> <output_file>"
    exit 1
fi

# Assign command-line arguments to variables
enhancers_bed=$1
gff3_file=$2
output_file=$3

# Add header line to the output file
echo -e "chr\tstart\tend\tid\ttype\tgene_intersection\tbp_to_closest_gene\tlength" > "$output_file"

# Append the processed data to the output file
bedtools closest -a <(sort -k1,1 -k2,2n "$enhancers_bed") \
    -b <(awk '$3 == "gene"' "$gff3_file" | sort -k1,1 -k4,4n) -d -t first | \
awk 'BEGIN{OFS="\t"} {
    dist = $NF; 
    len = $3 - $2;
    label = (dist == 0 ? "intronic" : "intergenic"); 
    print $1, $2, $3, $4, $5, label, dist, len
}' >> "$output_file"