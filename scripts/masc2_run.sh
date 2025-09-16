#!/bin/bash
#SBATCH -n 1
#SBATCH --time=12:0:0
#SBATCH --mem=64G
#SBATCH --job-name=macs2

bam_file=$1
out_dir=$2
name=$3
macs2 callpeak -t $bam_file -f BAM -n $name --outdir $out_dir -q 0.05
macs2 --version
