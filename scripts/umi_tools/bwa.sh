#!/bin/bash
cores=$1
ref_genome=$2
fastq=$3
output_bam=$4

bwa mem -t "$cores" "$ref_genome" "$fastq" | samtools view -Sb - > "$output_bam"