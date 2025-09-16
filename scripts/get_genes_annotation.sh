#!/bin/bash
#SBATCH -n 1
#SBATCH --time=1:0:0
#SBATCH --mem=8G

# download from: https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/
out_dir=$1
out_merged_file=$2
wget https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz -P ${out_dir}
gunzip ${out_dir}Homo_sapiens.GRCh38.113.gff3.gz

wget https://ftp.ensembl.org/pub/release-113/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v113.gff3.gz -P ${out_dir}
gunzip ${out_dir}Homo_sapiens.GRCh38.regulatory_features.v113.gff3.gz

cat ${out_dir}Homo_sapiens.GRCh38.113.gff3 >$out_merged_file
tail -n +2 ${out_dir}Homo_sapiens.GRCh38.regulatory_features.v113.gff3 >> $out_merged_file