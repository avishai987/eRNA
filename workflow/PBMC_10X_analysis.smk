module PBMC_10X_RNA:
    snakefile: "PBMC_10X_RNA.smk"


use rule * from PBMC_10X_RNA as PBMC_10X_RNA_*

module PBMC_10X_ATAC:
    snakefile: "PBMC_10X_ATAC.smk"
use rule * from PBMC_10X_ATAC as PBMC_10X_ATAC_*