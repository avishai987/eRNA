#!/bin/bash

#SBATCH -n 1
#SBATCH --time=48:0:0
#SBATCH --mem=2G
#SBATCH -J fasterq-dump
#SBATCH --output=GSE126074_SNARE_seq/logs/download_fastq_%j.log
. $lmod

module load sratoolkit
output_folder=$1
temp_data_folder=$2
SRR_number=$3
fasterq-dump $SRR_number --progress --outdir $output_folder --temp $temp_data_folder 

fasterq-dump --version