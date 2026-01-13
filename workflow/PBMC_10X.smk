
# ---------------------------------------------------------------------------- #
#                               10X_PBMC_GEX_ATAC                              #
# ---------------------------------------------------------------------------- #
#data description link: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k


# ------------------------------RNA Process------------------------------------- #
rule download_gex_matrix:
    params:
        dir = "Analysis/10X_PBMC/01_raw_data/gex_matrix/"
    output:
        log = "Analysis/10X_PBMC/01_raw_data/gex_matrix/log.txt",
        barcodes = "Analysis/10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
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
        log = "Analysis/10X_PBMC/01_raw_data/download_fastq_log.txt",
        fastq = "Analysis/10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k_fastqs.tar"
    conda:
        "eRNA"
    params:
        dir = "Analysis/10X_PBMC/01_raw_data/"
    shell: 
        '''
        mkdir -p {params.dir};
        {input.script} \
        'https://cg.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_fastqs.tar' \
        {params.dir}
        > {output.log} 2>&1;
        '''

        
# untar pbmc_granulocyte_sorted_10k_fastqs.tar
rule untar_fastq:
    input:
        fastq_tar = rules.download_fastq.output.fastq
    output:
        log = "Analysis/10X_PBMC/01_raw_data/untar_fastq_log.txt",
        fastq_rna_files = expand(
            "Analysis/10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_{lane}_{read}_001.fastq.gz",
            lane=["L003", "L004"],
            read=["R1", "R2"]
        )
    params:
        dir = "Analysis/10X_PBMC/01_raw_data/"
    shell: 
        '''
        tar -xvf {input.fastq_tar} --directory {params.dir} \
        > {output.log} 2>&1;
        '''

rule concatenation:
    input:
        R1_L3 = "Analysis/10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L003_R1_001.fastq.gz",
        R1_L4 = "Analysis/10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L004_R1_001.fastq.gz",
        R2_L3 = "Analysis/10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L003_R2_001.fastq.gz",
        R2_L4 = "Analysis/10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/gex/pbmc_granulocyte_sorted_10k_S1_L004_R2_001.fastq.gz"
    output:
        R1 = "Analysis/10X_PBMC/01_raw_data/concat_fastq/pbmc_granulocyte_sorted_10k_R1.fastq.gz",
        R2 = "Analysis/10X_PBMC/01_raw_data/concat_fastq/pbmc_granulocyte_sorted_10k_R2.fastq.gz",
    params:
        dir = "Analysis/10X_PBMC/01_raw_data/concat_fastq/"
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
        whitelist = rules.download_gex_matrix.output.barcodes
    output:
        fastq_R2_extracted = "Analysis/10X_PBMC/02_extract/pbmc_granulocyte_sorted_10k_R2_extracted.fastq.gz"
    params:
        dir = "Analysis/10X_PBMC/02_extract/",
        bc_pattern = "CCCCCCCCCCCCCCCCNNNNNNNNNNNN" # from https://assets.ctfassets.net/an68im79xiti/3o5BVs85etmhscbOF4XpP3/e4ad4cacf0a5b163a0c0e1f37d7b5cf7/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf
    conda:
        "eRNA_umi_tools"
    log:
        "Analysis/10X_PBMC/02_extract/extract.log"
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
        bam = "Analysis/10X_PBMC/03_align/pbmc_granulocyte_sorted_10k.bam"
    params:
        dir = "Analysis/10X_PBMC/03_align/",
        cores = 8,
        ref_genome = config["hg_38_bwa_index"]
    conda:
        "eRNA_bwa"
    log:
        "Analysis/10X_PBMC/03_align/bwa.log"
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
        counts = "Analysis/10X_PBMC/04_count/pbmc_granulocyte_sorted_10k_counts.txt",
        bam = "Analysis/10X_PBMC/04_count/pbmc_granulocyte_sorted_10k.bam.no_XS.bam.featureCounts.bam",
    params:
        dir = "Analysis/10X_PBMC/04_count/",
        feature_type = "regulatory_region",
        threads = 8
    conda:
        "eRNA_subread"
    log:
        "Analysis/10X_PBMC/04_count/featureCounts.log"
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
        sorted_tagged_bam = "Analysis/10X_PBMC/04_count/pbmc_granulocyte_sorted_10k_tagged_sorted.bam",
        bai_tagged_bam = "Analysis/10X_PBMC/04_count/pbmc_granulocyte_sorted_10k_tagged_sorted.bam.bai"
    params:
        dir = "Analysis/10X_PBMC/04_count/",
        threads = 8
    conda:
        "eRNA_samtools"
    log:
        "Analysis/10X_PBMC/04_count/sort_index.log"
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
        count_matrix = "Analysis/10X_PBMC/04_count/pbmc_granulocyte_sorted_10k_counts_per_cell.txt",
    params:
        dir = "Analysis/10X_PBMC/04_count/"
    conda:
        "eRNA_umi_tools"
    log:
        "Analysis/10X_PBMC/04_count/count_per_cell.log"
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
        whitelist = "Analysis/10X_PBMC/ATAC/02_extract/whitelist.tsv"
    params:
        dir = "Analysis/10X_PBMC/ATAC/02_extract/"
    conda:
        "eRNA_samtools"
    log:
        "Analysis/10X_PBMC/ATAC/02_extract/whitelist.log"
    shell:
        """
        mkdir -p {params.dir};
        # unzip whitelist and remove suffix after '-' to match cell barcodes in read1
        gzip -d -c {input.barcodes} | cut -f 1 -d '-' > {output.whitelist} 2> {log};
        """

rule download_fragments:
    output:
        fragments = "Analysis/10X_PBMC/01_raw_data/ATAC/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz"
    params:
        url = "https://cf.10xgenomics.com/samples/cell-arc/1.0.0/pbmc_granulocyte_sorted_10k/pbmc_granulocyte_sorted_10k_atac_fragments.tsv.gz",
        dir = "Analysis/10X_PBMC/01_raw_data/ATAC/"
    log:
        "Analysis/10X_PBMC/01_raw_data/ATAC/download_fragments.log"
    shell:
        '''
        wget  -P {params.dir} {params.url} -o {log}
        '''
rule count_fragments:
    input:
        fragments = rules.download_fragments.output.fragments,
        unique_enhancers = config["enhancers_to_count"]["bed"],
        count_by_fragments_script = "scripts/count_by_fragments.r",
        whitelist = "Analysis/10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        counts = "Analysis/10X_PBMC/ATAC/01_count_fragments/pbmc_granulocyte_sorted_10k_enhancer_counts.mtx",
        features = "Analysis/10X_PBMC/ATAC/01_count_fragments/pbmc_granulocyte_sorted_10k_features.tsv",
        barcodes = "Analysis/10X_PBMC/ATAC/01_count_fragments/pbmc_granulocyte_sorted_10k_barcodes.tsv"
    params:
        project_name = "pbmc_granulocyte_sorted_10k",
        out_dir = "Analysis/10X_PBMC/ATAC/01_count_fragments/",
        threads = 12
    conda:
        "erna_r"
    log:
        "Analysis/10X_PBMC/ATAC/01_count_fragments/count_fragments.log"
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

# ------------------------------PBMC enhancers------------------------------------- #

rule download_pbmc_enhancers:
    input:
        enhancers_file = "Analysis/10X_PBMC/cell_type_enhancers/enhancer_links.txt",
        combine_enhancers_script = "scripts/combine_enhancers_per_cellType.sh"
    output:
        combined_enhancers_loci = "Analysis/10X_PBMC/cell_type_enhancers/combined_enhancers_loci.bed"
    params:
        dir = "Analysis/10X_PBMC/cell_type_enhancers/"
    log:
        "Analysis/10X_PBMC/cell_type_enhancers/download_pbmc_enhancers.log"
    shell:
        '''
        mkdir -p {params.dir};
        wget  -i {input.enhancers_file} -P {params.dir}  > {log} 2>&1;
        {input.combine_enhancers_script} {output.combined_enhancers_loci} {params.dir} >> {log} 2>&1;
        '''
        
rule intersect_with_unique_enhancers:
    input:
        pbmc_enhancers = rules.download_pbmc_enhancers.output.combined_enhancers_loci,
        unique_enhancers = config["enhancers_to_count"]["bed"]
    output:
        intersected_enhancers = "Analysis/10X_PBMC/cell_type_enhancers/PBMC_enhancers_id.tsv"
    conda:
        "eRNA_bedtools"
    params:
        fraction = 0.5,
        dir = "Analysis/10X_PBMC/cell_type_enhancers/"
    shell:
        '''
        #remove 'chr' prefix
        sed -i 's/^chr//g' {input.pbmc_enhancers};
        bedtools intersect -a {input.unique_enhancers} -b {input.pbmc_enhancers} -f {params.fraction}  -wb | cut -f 4,9 > {output.intersected_enhancers};
        '''

# ------------------------------eRNA data analysis------------------------------------- #

#render scripts/PBMC_10K/gex_clustering.ipynb
rule gex_clustering:
    input:
        script = "scripts/PBMC_10K/gex_clustering.ipynb",
        mtx_file = "Analysis/10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/matrix.mtx.gz",
        features_file = "Analysis/10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/features.tsv.gz",
        barcodes_file = "Analysis/10X_PBMC/01_raw_data/gex_matrix/filtered_feature_bc_matrix/barcodes.tsv.gz"
    output:
        nb_out = "Analysis/10X_PBMC/gex_clustering/gex_clustering.ipynb",
        report = "Analysis/10X_PBMC/gex_clustering/gex_clustering.html",
        idents_out = "Analysis/10X_PBMC/gex_clustering/PBMC_RNA_idents.RDS",
        lib_size_out = "Analysis/10X_PBMC/gex_clustering/PBMC_RNA_lib_size.RDS"
    conda:
        "eRNA_jupyter"
    params:
        dir = "Analysis/10X_PBMC/gex_clustering/"
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
        rna_enhancers_counts_path = "Analysis/10X_PBMC/04_count/pbmc_granulocyte_sorted_10k_counts_per_cell.txt",
        enhancers_metadata = "Analysis/enhancers/ensembl/ensembl_enhancers_metadata.txt"
    output:
        nb_out = "Analysis/10X_PBMC/05_erna_preprocess/erna_preprocess.ipynb",
        report = "Analysis/10X_PBMC/05_erna_preprocess/erna_preprocess.html",
        filtered_erna = "Analysis/10X_PBMC/05_erna_preprocess/filtered_erna_pbmc_granulocyte_sorted_10k.rds"
    conda:
        "eRNA_jupyter"
    params:
        dir = "Analysis/10X_PBMC/05_erna_preprocess/"
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
        atac_metadata = "Analysis/10X_PBMC/03_ATAC/pbmc_granulocyte_sorted_10k_atac_metadata.txt",
        atac_enhancer_counts = "Analysis/10X_PBMC/03_ATAC/pbmc_granulocyte_sorted_10k_atac_counts_per_enhancer.txt"
    output:
        nb_out = "Analysis/10X_PBMC/06_analyze_by_ATAC/analyze_by_ATAC.ipynb",
        report = "Analysis/10X_PBMC/06_analyze_by_ATAC/analyze_by_ATAC.html"
    conda:
        "eRNA_jupyter"
    params:
        dir = "Analysis/10X_PBMC/06_analyze_by_ATAC/"
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
