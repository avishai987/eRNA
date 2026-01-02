#!/usr/bin/env bash

WHITELIST="$1"
FASTQ="$2"
out_file="$3"
# Filter single-end FASTQ reads by a whitelist of barcodes. Extracts barcode from the header (@BC:...),
#keeps reads whose barcode appears in whitelist, and writes them compressed output file

# explanation of the code:
#  paste - - - -        # Combine each read's 4 lines into one line
#awk -v WL="$whitelist" '                    # Pass whitelist file path to awk
#NR==FNR { wl[$1]; next }                # Step 1: load whitelist barcodes into hash "wl"
#            split($1, a, ":");                  # Split header by ":"
#            bc = substr(a[1], 2);               # Remove "@" to extract barcode
#            keep = (bc in wl);                  # Check if barcode is in whitelist
#            if (keep) print $0;                 # Print the whole read (4 lines in one)
#        }' "$whitelist" - | \
#    tr "\t" "\n" | gzip > "$out_file"           # Convert back to FASTQ 4-line format and compress


zcat "$FASTQ" | \
paste - - - - | \
awk -v WL="$WHITELIST" '
    NR==FNR { wl[$1]; next }
    {
        split($1, a, ":");
        bc = substr(a[1], 2);
        echo "sdf";
        keep = (bc in wl);
        if (keep) print $0;
    }' "$WHITELIST" - | \
tr "\t" "\n" | gzip > "$out_file"

