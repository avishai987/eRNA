#!/bin/bash

#SBATCH -n 1
#SBATCH --time=12:0:0
#SBATCH --mem=15G
#SBATCH -J tsv_to_h5

. $lmod
py_script=$1
path_to_rna_tsv=$2
path_to_atac_tsv=$3
output_file=$4
bash
export PATH="$conda_path:$PATH"
eval "$(conda shell.bash hook)"
conda activate h5py
python3 $py_script $path_to_rna_tsv $path_to_atac_tsv $output_file