#!/bin/bash

#SBATCH -n 4
#SBATCH --time=4:0:0
#SBATCH --mem=64G
#SBATCH -J run_RF
#SBATCH --mail-user=avishai.wizel@mail.huji.ac.il
#SBATCH --mail-type=END
conda_path=/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/bin/
lmod=/etc/profile.d/huji-lmod.sh
. $lmod
bash
export PATH="$conda_path:$PATH"
eval "$(conda shell.bash hook)"
conda activate tf_gpu_env
jupyter nbconvert --execute "/sci/labs/yotamd/lab_share/avishai.wizel/eRNA/10X_PBMC/scripts/05_random_forest.ipynb" --to 'html' --output "./05_random_forest.html"