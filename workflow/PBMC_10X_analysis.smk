module PBMC_10X_RNA:
    snakefile: "PBMC_10X_RNA.smk"


use rule * from PBMC_10X_RNA as PBMC_10X_RNA_*

module PBMC_10X_ATAC:
    snakefile: "PBMC_10X_ATAC.smk"
use rule * from PBMC_10X_ATAC as PBMC_10X_ATAC_*


module snare_star:
    snakefile: "snare.smk"

use rule * from snare_star

rule download_pbmc_enhancers:
    input:
        enhancers_file = "10X_PBMC/enhancers/enhancer_links.txt"
    output:
        B_cell_blood = "10X_PBMC/enhancers/B_cell_blood.bed",
        CD4 = "10X_PBMC/enhancers/CD4+.bed",
        CD8 = "10X_PBMC/enhancers/CD8+.bed",
        CD19 = "10X_PBMC/enhancers/CD19+.bed",
        CD20 = "10X_PBMC/enhancers/CD20+.bed",
        CD14 = "10X_PBMC/enhancers/CD14+.bed",
        CD14_monocyte = "10X_PBMC/enhancers/CD14+_monocyte.bed"
    params:
        dir = "10X_PBMC/enhancers/"
    log:
        "10X_PBMC/enhancers/download_pbmc_enhancers.log"
    shell:
        '''
        mkdir -p {params.dir};
        wget  -i {input.enhancers_file} -P {params.dir} -o {log};
        '''
        
rule intersect_with_unique_enhancers:
    input:
        pbmc_enhancers = "10X_PBMC/enhancers/{cell_type}.bed",
        unique_enhancers = rules.snare_07_unique_enhancers.output.uniqe_enhancers
    output:
        intersected_enhancers = "10X_PBMC/enhancers/enhancer_per_cell_type/{cell_type}_enhancers.tsv"
    conda:
        "eRNA_bedtools"
    params:
        fraction = 0.5,
        dir = "10X_PBMC/enhancers/enhancer_per_cell_type/"
    shell:
        '''
        #remove 'chr' prefix
        sed -i 's/^chr//g' {input.pbmc_enhancers};
        mkdir -p {params.dir};
        bedtools intersect -a {input.unique_enhancers} -b {input.pbmc_enhancers} -f {params.fraction} | cut -f 4 > {output.intersected_enhancers};
        '''

rule intersect_all_pbmc_enhancers:
    input:
        expand("10X_PBMC/enhancers/enhancer_per_cell_type/{cell_type}_enhancers.tsv", 
               cell_type=["B_cell_blood", "CD4+", "CD8+", "CD19+", "CD20+", "CD14+", "CD14+_monocyte"])