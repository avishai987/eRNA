
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

#render scripts/PBMC_10K/gex_clustering.ipynb
rule gex_clustering:
    input:
        script = "scripts/PBMC_10K/gex_clustering.ipynb",
        mtx_file = "10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/matrix.mtx.gz",
        features_file = "10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/features.tsv.gz",
        barcodes_file = "10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        nb_out = "10X_PBMC/gex_clustering/gex_clustering.ipynb",
        report = "10X_PBMC/gex_clustering/gex_clustering.html",
        idents_out = "10X_PBMC/gex_clustering/PBMC_RNA_idents.RDS",
        lib_size_out = "10X_PBMC/gex_clustering/PBMC_RNA_lib_size.RDS"
    conda:
        "eRNA_jupyter"
    params:
        dir = "10X_PBMC/gex_clustering/"
    shell:
        '''
        mkdir -p {params.dir};
        papermill {input.script} {output.nb_out} \
        -p mtx_file {input.mtx_file} \
        -p features_file {input.features_file} \
        -p barcodes_file {input.barcodes_file} \
        -p idents_out_path {output.idents_out} \
        -p lib_size_out_path {output.lib_size_out} \
        -k R --language R \
        && jupyter nbconvert --to html {output.nb_out} --output-dir {params.dir} 
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
# render scripts/PBMC_10K/analyze_by_ATAC.ipynb
rule analyze_by_ATAC:
    input:
        script = "scripts/PBMC_10K/analyze_by_ATAC.ipynb",
        filtered_erna = rules.erna_preprocess.output.filtered_erna,
        atac_metadata = "10X_PBMC/03_ATAC/pbmc_granulocyte_sorted_10k_atac_metadata.txt",
        atac_enhancer_counts = "10X_PBMC/03_ATAC/pbmc_granulocyte_sorted_10k_atac_counts_per_enhancer.txt"
    output:
        nb_out = "10X_PBMC/06_analyze_by_ATAC/analyze_by_ATAC.ipynb",
        report = "10X_PBMC/06_analyze_by_ATAC/analyze_by_ATAC.html"
    conda:
        "eRNA_jupyter"
    params:
        dir = "10X_PBMC/06_analyze_by_ATAC/"
    shell:
        '''
        mkdir -p {params.dir};
        papermill {input.script} {output.nb_out} \
        -p filtered_erna_path {input.filtered_erna} \
        -p atac_metadata_path {input.atac_metadata} \
        -p atac_enhancer_counts_path {input.atac_enhancer_counts} \
        -k R --language R \
        && jupyter nbconvert --to html {output.nb_out} --output-dir {params.dir} 
        '''
