

hg_38_bwa_index = config["hg_38_bwa_index"]
rule get_whitelist:
    input:
        barcodes = config["PBMC_10X"]["cells_whitelist"]
    output:
        whitelist = "10X_PBMC/ATAC/02_extract/whitelist.tsv"
    params:
        dir = "10X_PBMC/ATAC/02_extract/"
    conda:
        "eRNA_samtools"
    log:
        "10X_PBMC/ATAC/02_extract/whitelist.log"
    shell:
        """
        mkdir -p {params.dir};
        # unzip whitelist and remove suffix after '-' to match cell barcodes in read1
        gzip -d -c {input.barcodes} | cut -f 1 -d '-' > {output.whitelist} 2> {log};
        """

rule download_fragments:
    output:
        fragments = "10X_PBMC/01_raw_data/ATAC/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
    params:
        url = "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        dir = "10X_PBMC/01_raw_data/ATAC/"
    log:
        "10X_PBMC/01_raw_data/ATAC/download_fragments.log"
    shell:
        '''
        wget  -P {params.dir} {params.url} -o {log}
        '''
rule count_fragments:
    input:
        fragments = rules.download_fragments.output.fragments,
        unique_enhancers = config["enhancers_to_count"]["encode_bed"],
        count_by_fragments_script = "scripts/count_by_fragments.r",
        whitelist = "10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        counts = config["PBMC_10X"]["atac_matrix_counts"],
        features = config["PBMC_10X"]["atac_matrix_features"],
        barcodes = config["PBMC_10X"]["atac_matrix_barcode"]
    params:
        project_name = "pbmc_granulocyte_sorted_10k",
        out_dir = "10X_PBMC/ATAC/01_count_fragments/",
        threads = 12
    conda:
        "erna_r"
    log:
        "10X_PBMC/ATAC/01_count_fragments/count_fragments.log"
    shell:
        """
        Rscript {input.count_by_fragments_script} \
        {input.fragments} \
        {params.project_name} \
        {input.unique_enhancers} \
        {input.whitelist} \
        {params.out_dir} \
        {output.counts} \
        {params.threads}  > {log} 2>&1
        """ 