#!/bin/bash
#SBATCH -n 1
#SBATCH --time=3:0:0
#SBATCH --mem=30G
#SBATCH -J bedtools
. $lmod
module load samtools
bam_file=$1
filtered_barcode_list=$2
bam_out=$3


samtools view -D CB:$filtered_barcode_list $bam_file -o  $bam_out

module load bedtools2
samtools index $bam_out