#!/bin/bash

#SBATCH -n 1
#SBATCH --time=12:0:0
#SBATCH --mem=20G
#SBATCH -J BABEL.py

. $lmod
py_script=$1
path_to_rna_atac_h5=$2
our_dir=$3
bash
export PATH="$conda_path:$PATH"
eval "$(conda shell.bash hook)"
conda activate h5py
python3 $py_script --snareseq $path_to_rna_atac_h5 --outdir $our_dir

