# from https://assets.ctfassets.net/an68im79xiti/3o5BVs85etmhscbOF4XpP3/e4ad4cacf0a5b163a0c0e1f37d7b5cf7/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevA.pdf
# R1: 50 bp dna read 
# R2: 16 bp cell barcode , no UMI
# R3: 49 bp dna read 
rule extract_lane:
    input:
        fastq_R1 = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/atac/pbmc_granulocyte_sorted_10k_S16_{lane}_R1_001.fastq.gz",
        fastq_R2 = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/atac/pbmc_granulocyte_sorted_10k_S16_{lane}_R2_001.fastq.gz",
        fastq_R3 = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/atac/pbmc_granulocyte_sorted_10k_S16_{lane}_R3_001.fastq.gz",
        whitelist = "10X_PBMC/ATAC/02_extract/whitelist.tsv",
        filter_barcodes_script = "scripts/filter_barcodes.sh"
    output:
        fastq_R1_extracted = "10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_R1_001.barcoded.fastq.gz",
        fastq_R3_extracted = "10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_R3_001.barcoded.fastq.gz",
    params:
        dir = "10X_PBMC/ATAC/02_extract/",
        fastq_dir = "10X_PBMC/01_raw_data/pbmc_granulocyte_sorted_10k/atac/",
        bc_length = 16
    conda:
        "eRNA_sinto"
    log:
        "10X_PBMC/ATAC/02_extract/extract_{lane}.log"
    shell:
        """
        # extract cell barcodes from R2 and add to R1 and R3:
        mkdir -p {params.dir};
        sinto barcode -b {params.bc_length}  --barcode_fastq {input.fastq_R2} \
        --read1 {input.fastq_R1}  \
        --read2 {input.fastq_R3}  \
        --whitelist {input.whitelist}  > {log} 2>&1;
        # move output files to desired location
        mv {params.fastq_dir}pbmc_granulocyte_sorted_10k_S16_{wildcards.lane}_R1_001.barcoded.fastq.gz {params.dir};
        mv {params.fastq_dir}pbmc_granulocyte_sorted_10k_S16_{wildcards.lane}_R3_001.barcoded.fastq.gz {params.dir};
        """

LANES = ["L001", "L002", "L003","L004"]

rule extract_all_lanes:
    input:
        expand("10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_R1_001.barcoded.fastq.gz", 
               lane=LANES),
        expand("10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_R3_001.barcoded.fastq.gz", 
               lane=LANES)
# rule to filter fastq with scripts/filter_barcodes.sh
rule filter_barcodes:
    input:
        whitelist = "10X_PBMC/ATAC/02_extract/whitelist.tsv",
        fastq_in = "10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_{read}_001.barcoded.fastq.gz"
    output:
        fastq_out = "10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_{read}_001.barcoded_filtered.fastq.gz"
    params:
        script = "scripts/filter_barcodes.sh"
    shell:
        """
        {params.script} {input.whitelist} {input.fastq_in} {output.fastq_out}
        """

rule filter_all_reads:
    input:
        expand("10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_{read}_001.barcoded_filtered.fastq.gz", 
               lane=LANES, read=["R1", "R3"])
               
rule align_lane:
    input:
        fastq_R1_extracted = "10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_R1_001.barcoded.fastq.gz",
        fastq_R3_extracted = "10X_PBMC/ATAC/02_extract/pbmc_granulocyte_sorted_10k_S16_{lane}_R3_001.barcoded.fastq.gz",
    output:
        bam = "10X_PBMC/ATAC/03_align/pbmc_granulocyte_sorted_10k_S16_{lane}_aligned.bam"
    params:
        dir = "10X_PBMC/ATAC/03_align/",
        ref_genome = config["hg_38_bwa_index"]

    conda:
        "eRNA_bwa"
    log:
        "10X_PBMC/ATAC/03_align/align_{lane}.log"
    shell:
        """
        mkdir -p {params.dir};
        bwa mem -t 8 {params.ref_genome} {input.fastq_R1_extracted} {input.fastq_R3_extracted} 2> {log} | \
        samtools view -bS - > {output.bam}
        """

rule align_all_lanes:
    input:
        expand("10X_PBMC/ATAC/03_align/pbmc_granulocyte_sorted_10k_S16_{lane}_aligned.bam", 
               lane=LANES)
