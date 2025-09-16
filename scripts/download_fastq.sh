#!/bin/bash

#SBATCH -n 1
#SBATCH --time=48:0:0
#SBATCH --mem=12G
#SBATCH -J fasterq-dump
#SBATCH --output=GSE126074_SNARE_seq/logs/download_fastq_%j.log
#SBATCH --mail-type=END
#SBATCH --mail-user=avishai.wizel@mail.huji.ac.il

output_folder=$1
temp_data_folder=$2
SRR_number=$3 

fasterq-dump $SRR_number --progress --outdir $output_folder --temp $temp_data_folder 
# gzip the files
gzip $output_folder/${SRR_number}_1.fastq
# if paired-end
if [ -f $output_folder/${SRR_number}_2.fastq ]; then
  gzip $output_folder/${SRR_number}_2.fastq
fi

fasterq-dump --version