#!/bin/bash

#SBATCH -n 8
#SBATCH --time=12:0:0
#SBATCH --mem=30G
#SBATCH -J star
#SBATCH --output=GSE126074_SNARE_seq/logs/alignment_%j.log


genomeDir=$1
fastq1=$2
fastq2=$3
output_prefix=$4

STAR --genomeDir $genomeDir --readFilesIn $fastq2 $fastq1 --readFilesCommand zcat --soloCBwhitelist \
  None --soloCellFilter EmptyDrops_CR 3000 0.99 10 45000 90000 500 0.01 20000 0.001 10000 --outFileNamePrefix $output_prefix \
  --runThreadN 8 --clip3pAdapterSeq AAAAAAAAAAAA --outFilterMatchNmin 50 --soloFeatures Gene --outSAMattributes NH HI nM AS CR UR CB UB GX GN \
  --soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 9 --soloCBmatchWLtype 1MM --outSAMtype BAM SortedByCoordinate \
  --soloCellReadStats Standard --soloBarcodeReadLength 0

STAR --version