#!/bin/bash
genes_gtf=$1
bam_file=$2
output_count=$3
threads=$4
feature_type=$5
featureCounts -a "$genes_gtf" -o "$output_count" -R BAM "$bam_file" -T "$threads" -t "$feature_type"
