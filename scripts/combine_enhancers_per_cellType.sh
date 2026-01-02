#!/bin/bash

# Define the final output filename
output_file=$1
workdir=$2
# Initialize/Empty the output file if it already exists
> "$output_file"

# Loop through all BED files in the current directory
for file in "$workdir"/*.bed; do
    # Skip the output file itself if it ends in .bed to avoid infinite loops
    if [[ "$file" == "$output_file" ]]; then continue; fi

    # 1. Extract the cell type name from the filename (e.g., 'T-cells.bed' -> 'T-cells')
    cell_type=$(basename "$file" .bed)
    
    echo "Processing cell type: $cell_type"
    
    # 2. Use AWK to reformat the columns into a BED6 standard:
    # $1, $2, $3: Chromosome, Start, End (Original coordinates)
    # ct (cell_type): Moved to Column 4 (Name field)
    # $4 (enrichment score): Moved to Column 5 (Score field)
    # ".": Added to Column 6 (Strand field as a placeholder)
    awk -v ct="$cell_type" 'BEGIN {OFS="\t"} {print $1, $2, $3, ct, $4, "."}' "$file" >> "$output_file"
done

# 3. Sort the combined file by genomic coordinates (Chr then Start position)
# This is mandatory for downstream tools like bedtools and tabix
echo "Sorting the combined file..."
sort -k1,1 -k2,2n "$output_file" > "$workdir/combined_cells.sorted.bed"

# 4. Optional: Clean up the unsorted file
rm "$output_file"
#rename the sorted file to the original output filename
mv "$workdir/combined_cells.sorted.bed" "$output_file"
echo "Done! Final file: $output_file"