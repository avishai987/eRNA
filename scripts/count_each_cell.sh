#!/bin/bash
source ~/.bashrc
module load py-htseq

bam_file=$1
gtf_file=$2
output=$3
htseq-count-barcodes $bam_file $gtf_file -t regulatory_region --counts_output $output --stranded no