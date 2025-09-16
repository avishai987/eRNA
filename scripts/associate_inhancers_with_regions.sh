#! /bin/bash
# Associate enhancers with regions
regions_file=$1
enhancers_file=$2
output_file=$3

bedtools intersect -a $regions_file -b $enhancers_file -wao > $output_file
