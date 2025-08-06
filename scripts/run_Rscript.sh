#!/bin/bash
#SBATCH -n 1
#SBATCH -J Rscript
. $lmod
module load R4/4.4.1
Rscript -e "rmarkdown::render(
    input = '$1',
    output_format = 'html_document',
    output_file = '$2',
    knit_root_dir = getwd(),
    output_dir = dirname('$2')
  )
"