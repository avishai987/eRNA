#!/bin/zsh

#SBATCH -n 1
#SBATCH --time=12:0:0
#SBATCH --mem=16G
#SBATCH -J bedtools
atac_seq_file=$1
annotations_bed_file=$2
out_file=$3
#make bed from atac loci 
atac_bed=$(cat "$atac_seq_file" | awk 'NR > 1 {print $1}' | awk -F '[:-]' '{OFS="\t"; print $1, $2, $3}')
#intersect with enhacers annotations
echo $atac_bed | bedtools intersect -wa -a stdin -b $annotations_bed_file | awk '{printf "%s:%s-%s\n", $1, $2, $3}' > $out_file



