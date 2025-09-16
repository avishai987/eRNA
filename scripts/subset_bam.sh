#!/bin/bash
#SBATCH -n 1
#SBATCH --time=3:0:0
#SBATCH --mem=30G
#SBATCH -J bedtools

bam_file=$1
filtered_barcode_list=$2
bam_out=$3

echo $bam_file
echo $filtered_barcode_list
echo $bam_out
samtools view -D CB:$filtered_barcode_list $bam_file -o  $bam_out

samtools index $bam_out