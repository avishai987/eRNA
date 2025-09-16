#!/bin/bash -i
# Snakemake jobscript template

# load conda functions (אם use-conda:true)
source /sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/etc/profile.d/conda.sh

echo "Job started on $(date)"

# run the command
{exec_job}