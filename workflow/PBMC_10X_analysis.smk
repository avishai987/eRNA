
rule download_pbmc_enhancers:
    input:
        enhancers_file = "10X_PBMC/cell_type_enhancers/enhancer_links.txt",
        combine_enhancers_script = "scripts/combine_enhancers_per_cellType.sh"
    output:
        combined_enhancers_loci = "10X_PBMC/cell_type_enhancers/combined_enhancers_loci.bed"
    params:
        dir = "10X_PBMC/cell_type_enhancers/"
    log:
        "10X_PBMC/cell_type_enhancers/download_pbmc_enhancers.log"
    shell:
        '''
        mkdir -p {params.dir};
        wget  -i {input.enhancers_file} -P {params.dir}  > {log} 2>&1;
        {input.combine_enhancers_script} {output.combined_enhancers_loci} {params.dir} >> {log} 2>&1;
        '''
        
rule intersect_with_unique_enhancers:
    input:
        pbmc_enhancers = rules.download_pbmc_enhancers.output.combined_enhancers_loci,
        unique_enhancers = config["enhancers_to_count"]["encode_bed"]
    output:
        intersected_enhancers = "10X_PBMC/cell_type_enhancers/PBMC_enhancers_id.tsv"
    conda:
        "eRNA_bedtools"
    params:
        fraction = 0.5,
        dir = "10X_PBMC/cell_type_enhancers/"
    shell:
        '''
        #remove 'chr' prefix
        sed -i 's/^chr//g' {input.pbmc_enhancers};
        bedtools intersect -a {input.unique_enhancers} -b {input.pbmc_enhancers} -f {params.fraction}  -wb | cut -f 4,9 > {output.intersected_enhancers};
        '''

#render scripts/PBMC_10K/erna_preprocess.ipynb
rule erna_preprocess:
    input:
        script = "scripts/PBMC_10K/erna_preprocess.ipynb",
        rna_enhancers_counts_path = "10X_PBMC/04_count/pbmc_granulocyte_sorted_10k_counts_per_cell.txt",
        enhancers_metadata = "Analysis/enhancers/ensembl/ensembl_enhancers_metadata.txt"
    output:
        nb_out = "10X_PBMC/05_erna_preprocess/erna_preprocess.ipynb",
        report = "10X_PBMC/05_erna_preprocess/erna_preprocess.html",
        filtered_erna = "10X_PBMC/05_erna_preprocess/filtered_erna_pbmc_granulocyte_sorted_10k.rds"
    conda:
        "eRNA_jupyter"
    params:
        dir = "10X_PBMC/05_erna_preprocess/"
    shell:
        '''
        mkdir -p {params.dir};
        papermill {input.script} {output.nb_out} \
        -p min_cells 10 -p min_enhancer_counts 10 -p rna_enhancers_counts_path {input.rna_enhancers_counts_path}\
        -k R --language R \
        && jupyter nbconvert --to html {output.nb_out} --output-dir {params.dir} 
        '''

