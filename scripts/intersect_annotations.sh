#!/bin/bash
#SBATCH -n 1
#SBATCH --time=3:0:0
#SBATCH --mem=8G
#SBATCH -J bedtools
. $lmod
module load bedtools2
bed_1=$1
bed_2=$2
bed_3=$3
bed_4=$4
out_file=$5
bedtools intersect -a $bed_1 -b $bed_2 $bed_3 $bed_4 > $out_file