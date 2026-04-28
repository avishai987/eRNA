# ---------------------------------------------------------------------------- #
#                               LN_10X_GEX_ATAC                              #
# ---------------------------------------------------------------------------- #
#data description link: https://www.10xgenomics.com/datasets/fresh-frozen-lymph-node-with-b-cell-lymphoma-14-k-sorted-nuclei-1-standard-2-0-0


# ------------------------------RNA Process------------------------------------- #
rule download_gex_matrix:
    params:
        dir = "Analysis/10X_LN/01_raw_data/gex_matrix/"
    output:
        barcodes = "Analysis/10X_LN/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz",
        features = "Analysis/10X_LN/01_raw_data/gex_matrix/filtered_feature_bc_matrix/features.tsv.gz",
        matrix = "Analysis/10X_LN/01_raw_data/gex_matrix/filtered_feature_bc_matrix/matrix.mtx.gz"
    log:
        "Analysis/10X_LN/01_raw_data/gex_matrix/download_gex_matrix_log.txt"
    shell:
        '''
        mkdir -p {params.dir}
        wget https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.tar.gz \
         -P {params.dir} > {log} 2>&1;
         tar -xzvf {params.dir}/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.tar.gz -C {params.dir} >> {log} 2>&1;
         rm {params.dir}/lymph_node_lymphoma_14k_filtered_feature_bc_matrix.tar.gz
        '''

rule download_fastq:
    input:
        script = "scripts/wget_file.sh"
    output:
        log = "Analysis/10X_LN/01_raw_data/tar/download_fastq_log.txt",
        fastq = "Analysis/10X_LN/01_raw_data/tar/lymph_node_lymphoma_14k_fastqs.tar"
    conda:
        "eRNA"
    params:
        dir = "Analysis/10X_LN/01_raw_data/tar/"
    shell: 
        '''
        mkdir -p {params.dir};
        {input.script} \
        'https://s3-us-west-2.amazonaws.com/10x.files/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_fastqs.tar' \
        {params.dir}
        > {output.log} 2>&1;
        '''
rule untar_fastq:
    input:
        fastq_tar = rules.download_fastq.output.fastq
    output:
            fastq_rna_files = expand(
                "Analysis/10X_LN/01_raw_data/extracted_tar/fastqs/gex/lymph_node_lymphoma_14k_S1_L002_{read}_001.fastq.gz",
                read=["R1", "R2"]
            )
    params:
        dir = "Analysis/10X_LN/01_raw_data/extracted_tar/"
    log:
        "Analysis/10X_LN/01_raw_data/extracted_tar/untar_fastq_log.txt"
    shell: 
        '''
        tar -xvf {input.fastq_tar} --directory {params.dir} fastqs/gex/  \
        > {log} 2>&1;
        '''
rule extract:
    input:
        fastq_R1 = rules.untar_fastq.output.fastq_rna_files[0],
        fastq_R2 = rules.untar_fastq.output.fastq_rna_files[1],
        whitelist = rules.download_gex_matrix.output.barcodes
    output:
        fastq_R2_extracted = "Analysis/10X_LN/02_extract/lymph_node_lymphoma_14k_R2_extracted.fastq.gz"
    params:
        dir = "Analysis/10X_LN/02_extract/",
        bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN" # from https://assets.ctfassets.net/an68im79xiti/3o5BVs85etmhscbOF4XpP3/e4ad4cacf0a5b163a0c0e1f37d7b5cf7/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf
    conda:
        "eRNA_umi_tools"
    log:
        "Analysis/10X_LN/02_extract/extract.log"
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
        bam = "Analysis/10X_LN/03_align/lymph_node_lymphoma_14k.bam"
    params:
        dir = "Analysis/10X_LN/03_align/",
        cores = 8,
        ref_genome = config["hg_38_bwa_index"]
    conda:
        "eRNA_bwa"
    log:
        "Analysis/10X_LN/03_align/bwa.log"
    shell:
        """
        mkdir -p {params.dir};
        bwa mem -t {params.cores} {params.ref_genome} {input.fastq_R2_extracted} 2> {log} | samtools view -Sb - > {output.bam};
        """
rule assign_to_genes:
    input:
        bam = rules.align.output.bam,
        gtf = config["enhancers_to_count"]["gtf"]
    output:
        counts = "Analysis/10X_LN/04_count/lymph_node_lymphoma_14k_counts.txt",
        bam = "Analysis/10X_LN/04_count/lymph_node_lymphoma_14k.bam.no_XS.bam.featureCounts.bam",
    params:
        dir = "Analysis/10X_LN/04_count/",
        feature_type = "regulatory_region",
        threads = 8
    conda:
        "eRNA_subread"
    log:
        "Analysis/10X_LN/04_count/featureCounts.log"
    shell:
        """
        mkdir -p {params.dir};
        # remove XS tag from BAM (added by BWA), as it causes umi_tools to fail
        samtools view -h {input.bam} --remove-tag XS -b > {input.bam}.no_XS.bam;
        featureCounts -a {input.gtf} -o {output.counts} -R BAM {input.bam}.no_XS.bam -T {params.threads} -t {params.feature_type} > {log} 2>&1;
        rm -f {input.bam}.no_XS.bam {input.bam}.featureCounts.bam;
        """

rule sort_index:
    input:
        bam = rules.assign_to_genes.output.bam
    output:
        sorted_tagged_bam = "Analysis/10X_LN/04_count/lymph_node_lymphoma_14k_tagged_sorted.bam",
        bai_tagged_bam = "Analysis/10X_LN/04_count/lymph_node_lymphoma_14k_tagged_sorted.bam.bai"
    params:
        dir = "Analysis/10X_LN/04_count/",
        threads = 8
    conda:
        "eRNA_samtools"
    log:
        "Analysis/10X_LN/04_count/sort_index.log"
    shell:
        """
        mkdir -p {params.dir};
        samtools sort -@ {params.threads} -m 8G {input.bam} -o {output.sorted_tagged_bam} > {log} 2>&1;
        samtools index {output.sorted_tagged_bam} >> {log} 2>&1;
        """
rule count_per_cell:
    input:
        sorted_tagged_bam = rules.sort_index.output.sorted_tagged_bam
    output:
        count_matrix = "Analysis/10X_LN/04_count/lymph_node_lymphoma_14k_counts_per_cell.txt",
    params:
        dir = "Analysis/10X_LN/04_count/"
    conda:
        "eRNA_umi_tools"
    log:
        "Analysis/10X_LN/04_count/count_per_cell.log"
    shell:
        """
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS \
        --per-cell  -I {input.sorted_tagged_bam} -S {output.count_matrix} \
        --wide-format-cell-counts > {log} 2>&1;
        """


# ------------------------------ATAC Process------------------------------------- #


rule get_whitelist:
    input:
        barcodes = rules.download_gex_matrix.output.barcodes
    output:
        whitelist = "Analysis/10X_LN/ATAC/02_extract/whitelist.tsv"
    params:
        dir = "Analysis/10X_LN/ATAC/02_extract/"
    conda:
        "eRNA_bamtools"
    log:
        "Analysis/10X_LN/ATAC/02_extract/whitelist.log"
    shell:
        """
        mkdir -p {params.dir};
        # unzip whitelist and remove suffix after '-' to match cell barcodes in read1
        gzip -d -c {input.barcodes} | cut -f 1 -d '-' > {output.whitelist} 2> {log};
        """

rule download_fragments:
    output:
        fragments = "Analysis/10X_LN/ATAC/01_raw_data/lymph_node_lymphoma_14k_atac_fragments.tsv.gz"
    params:
        url = "https://cf.10xgenomics.com/samples/cell-arc/2.0.0/lymph_node_lymphoma_14k/lymph_node_lymphoma_14k_atac_fragments.tsv.gz",
        dir = "Analysis/10X_LN/ATAC/01_raw_data/"
    log:
        "Analysis/10X_LN/ATAC/01_raw_data/download_fragments.log"
    shell:
        '''
        wget  -P {params.dir} {params.url} -o {log}
        '''
rule count_fragments:
    input:
        fragments = rules.download_fragments.output.fragments,
        unique_enhancers = config["enhancers_to_count"]["bed"],
        count_by_fragments_script = "scripts/count_by_fragments.r",
        whitelist = "Analysis/10X_LN/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        counts = "Analysis/10X_LN/ATAC/03_count_fragments/LN_sorted_10k_enhancer_counts.mtx",
        features = "Analysis/10X_LN/ATAC/03_count_fragments/LN_sorted_10k_features.tsv",
        barcodes = "Analysis/10X_LN/ATAC/03_count_fragments/LN_sorted_10k_barcodes.tsv"
    params:
        project_name = "LN_sorted_10k",
        out_dir = "Analysis/10X_LN/ATAC/03_count_fragments/",
        threads = 12
    conda:
        "erna_r"
    log:
        "Analysis/10X_LN/ATAC/03_count_fragments/count_fragments.log"
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

rule gex_clustering:
    input:
        script = "scripts/10X_LN/gex_clustering.ipynb",
        mtx_file = "Analysis/10X_LN/01_raw_data/gex_matrix/filtered_feature_bc_matrix/matrix.mtx.gz",
        features_file = "Analysis/10X_LN/01_raw_data/gex_matrix/filtered_feature_bc_matrix/features.tsv.gz",
        barcodes_file = "Analysis/10X_LN/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        nb_out = "Analysis/10X_LN/gex_clustering/gex_clustering.ipynb",
        report = "Analysis/10X_LN/gex_clustering/gex_clustering.html",
        idents_out = "Analysis/10X_LN/gex_clustering/LN_RNA_idents.RDS",
        lib_size_out = "Analysis/10X_LN/gex_clustering/LN_RNA_lib_size.RDS"
    conda:
        "eRNA_jupyter"
    params:
        dir = "Analysis/10X_LN/gex_clustering/"
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