
# ---------------------------------------------------------------------------- #
#                               10X_PBMC_GEX_ATAC                              #
# ---------------------------------------------------------------------------- #

rule PBMC_10X_01_download_atac:
    #data description link: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k
    shell:
        '''
        mkdir -p ./10X_PBMC/01_raw_data
        cd ./10X_PBMC/01_raw_data
        wget https://github.com/wukevin/babel/raw/refs/heads/main/data/10x/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
        '''
rule PBMC_10X_02_download_fastq:
    input:
        script = "scripts/wget_file.sh"
    output:
        log = "10X_PBMC/01_raw_data/download_fastq_log.txt",
        fastq = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k_fastqs.tar"
    conda:
        "eRNA"
    params:
        dir = "10X_PBMC/01_raw_data/"
    shell: 
        '''
        mkdir -p {params.dir};
        {input.script} \
        'https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_fastqs.tar' \
        {params.dir}
        > {output.log} 2>&1;
        '''
#check md5sum
rule PBMC_10X_03_check_md5sum:
    input:
        script = "scripts/md5sum.sh",
        fastq = rules.PBMC_10X_02_download_fastq.output.fastq
    output:
        log = "10X_PBMC/01_raw_data/check_md5sum_log.txt"
    params:
        dir = "10X_PBMC/01_raw_data/"
    shell: 
        '''
        mkdir -p {params.dir};
        {input.script} \
        {input.fastq} \ 
        > {output.log} 2>&1;
        '''
#untar pbmc_granulocyte_sorted_10k_fastqs.tar
rule PBMC_10X_04_untar_fastq:
    input:
        fastq = rules.PBMC_10X_02_download_fastq.output.fastq
    output:
        log = "10X_PBMC/01_raw_data/untar_fastq_log.txt"
    params:
        dir = "10X_PBMC/01_raw_data/"
    shell: 
        '''
        tar -xvf {input.fastq} --directory {params.dir} \
        > {output.log} 2>&1;
        '''