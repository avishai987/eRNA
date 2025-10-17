
# ---------------------------------------------------------------------------- #
#                               10X_PBMC_GEX_ATAC                              #
# ---------------------------------------------------------------------------- #
#data description link: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k

hg_38_bwa_index = "/sci/data/reference_data/Homo_sapiens/Ensembl/GRCh38/Sequence/BWAIndex/genome.fa"


rule download_gex_matrix:
    params:
        dir = "10X_PBMC/01_raw_data/gex_matrix/"
    output:
        "10X_PBMC/01_raw_data/gex_matrix/log.txt"
    shell:
        '''
        mkdir -p {params.dir}
        wget https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz \
         -P {params.dir} > {output.log} 2>&1;
         tar -xzvf {params.dir}/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz -C {params.dir} >> {output.log} 2>&1;
         rm {params.dir}/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.tar.gz
        '''

rule download_fastq:
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
        fastq = rules.download_fastq.output.fastq
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
rule untar_fastq:
    input:
        fastq = rules.download_fastq.output.fastq
    output:
        log = "10X_PBMC/01_raw_data/untar_fastq_log.txt"
    params:
        dir = "10X_PBMC/01_raw_data/"
    shell: 
        '''
        tar -xvf {input.fastq} --directory {params.dir} \
        > {output.log} 2>&1;
        '''

rule concatenation:
    input:
        R1_L3 = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L003_R1_001.fastq.gz",
        R1_L4 = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L004_R1_001.fastq.gz",
        R2_L3 = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L003_R2_001.fastq.gz",
        R2_L4 = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L004_R2_001.fastq.gz"
    output:
        R1 = "10X_PBMC/01_raw_data/concat_fastq/pbmc_granulocyte_sorted_10k_R1.fastq.gz",
        R2 = "10X_PBMC/01_raw_data/concat_fastq/pbmc_granulocyte_sorted_10k_R2.fastq.gz",
    params:
        dir = "10X_PBMC/01_raw_data/concat_fastq/"
    shell:
        '''
        mkdir -p {params.dir};
        cat {input.R1_L3} {input.R1_L4} > {output.R1} ;
        cat {input.R2_L3} {input.R2_L4} > {output.R2} ;
        '''

rule extract:
    input:
        fastq_R1 = rules.concatenation.output.R1,
        fastq_R2 = rules.concatenation.output.R2,
        whitelist = "10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        fastq_R2_extracted = "10X_PBMC/02_extract/pbmc_granulocyte_sorted_10k_R2_extracted.fastq.gz"
    params:
        dir = "10X_PBMC/02_extract/",
        bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN" # from https://assets.ctfassets.net/an68im79xiti/3o5BVs85etmhscbOF4XpP3/e4ad4cacf0a5b163a0c0e1f37d7b5cf7/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf
    conda:
        "eRNA_umi_tools"
    log:
        "10X_PBMC/02_extract/extract.log"
    shell:
        """
        mkdir -p {params.dir};
        # unzip whitelist and remove suffix after '-' to match cell barcodes in read1
        gzip -d -c {input.whitelist} | cut -f 1 -d '-' > {params.dir}/whitelist.tsv;
        pigz -d -c {input.fastq_R1} | umi_tools extract --bc-pattern={params.bc_pattern} \
        --read2-in={input.fastq_R2} \
        --read2-stdout \
        --log={log} \
        --whitelist {params.dir}/whitelist.tsv | pigz -c -p 8 > {output.fastq_R2_extracted}
        """

rule align:
    input:
        fastq_R2_extracted = rules.extract.output.fastq_R2_extracted
    output:
        bam = "10X_PBMC/03_align/pbmc_granulocyte_sorted_10k.bam"
    params:
        dir = "10X_PBMC/03_align/",
        cores = 8,
        ref_genome = hg_38_bwa_index
    conda:
        "eRNA_bwa"
    log:
        "10X_PBMC/03_align/bwa.log"
    shell:
        """
        mkdir -p {params.dir};
        bwa mem -t {params.cores} {params.ref_genome} {input.fastq_R2_extracted} 2> {log} | samtools view -Sb - > {output.bam};
        samtools sort -o {output.bam} {output.bam} >> {log} 2>&1;
        samtools index {output.bam} >> {log} 2>&1;
        """