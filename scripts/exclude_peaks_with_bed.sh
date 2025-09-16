#!/bin/bash

module load bedtools2

peaks=$1
genes=$2
output_file=$3
if [[ -z "$output_file" ]]; then
  echo "Error: output_file is empty" >&2
  exit 1
fi

bedtools intersect -a "${peaks}" -b "${genes}" -v > $output_file
bedtools --version
