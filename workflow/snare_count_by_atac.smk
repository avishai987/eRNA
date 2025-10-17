# Count open chromatin snare reads in enhancers


module snare_star:
    snakefile: "snare.smk"

use rule * from snare_star

hg_38_bwa_index = "/sci/data/reference_data/Homo_sapiens/Ensembl/GRCh38/Sequence/BWAIndex/genome.fa"

rule download_fastq:
    input:
        script = "scripts/download_fastq.sh"
    output:
        fastq1 = "GSE126074_SNARE_seq/count_by_atac/01_raw_data/SRR8528319_1.fastq.gz",
        fastq2 = "GSE126074_SNARE_seq/count_by_atac/01_raw_data/SRR8528319_2.fastq.gz",
        fastq3 = "GSE126074_SNARE_seq/count_by_atac/01_raw_data/SRR8528319_3.fastq.gz"
    conda:
        "eRNA_sra-tools"
    log:
        "GSE126074_SNARE_seq/count_by_atac/01_raw_data/download_fastq.log"
    params: 
        dir = "GSE126074_SNARE_seq/count_by_atac/01_raw_data/",
        srr = "SRR8528319"
    shell:
        '''
        mkdir -p {params.dir};
        {input.script} \
        {params.dir} \
        GSE126074_SNARE_seq \
        {params.srr} \
        > {log} 2>&1;
        '''

rule extract:
    input:
        script = "scripts/umi_tools/extract.sh",
        bc_fastq = rules.download_fastq.output.fastq1,
        fastq2 = rules.download_fastq.output.fastq2,
        fastq3 = rules.download_fastq.output.fastq3,
        whitelist = rules.snare_02_download_counts.output.barcodes
    output:
        fastq2_extracted = "GSE126074_SNARE_seq/count_by_atac/02_extract/SRR8528318_2_extracted.fastq.gz",
        fastq3_extracted = "GSE126074_SNARE_seq/count_by_atac/02_extract/SRR8528318_3_extracted.fastq.gz"
    params:
        dir = "GSE126074_SNARE_seq/count_by_atac/02_extract/",
        bc_pattern = "CCCCCCCCCCCCNNNNNNNNN"
    conda:
        "eRNA_umi_tools"
    log:
        "GSE126074_SNARE_seq/count_by_atac/02_extract/extract.log"
    shell:
        """
        mkdir -p {params.dir};
        {input.script} {params.bc_pattern} {input.bc_fastq} {input.fastq2} /dev/null  \
        {output.fastq2_extracted} {input.whitelist} > {log} 2>&1;
        {input.script} {params.bc_pattern} {input.bc_fastq} {input.fastq3} /dev/null  \
        {output.fastq3_extracted} {input.whitelist} >> {log} 2>&1;
        """
rule align:
    input:
        script = "scripts/umi_tools/bwa.sh",
        fastq2 = rules.extract.output.fastq2_extracted,
        fastq3 = rules.extract.output.fastq3_extracted
    output:
        bam = "GSE126074_SNARE_seq/count_by_atac/03_bwa/SRR8528319.bam"
    params:
        dir = "GSE126074_SNARE_seq/count_by_atac/03_bwa",
        threads = 4,
        ref_genome = hg_38_bwa_index
    conda:
        "eRNA_bwa"
    log:
        "GSE126074_SNARE_seq/count_by_atac/03_bwa/bwa.log"
    shell:
        """
        mkdir -p {params.dir};
        bwa mem -t {threads} {params.ref_genome} {input.fastq2} {input.fastq3} \
        | samtools view -b - > {output.bam}
        """ 