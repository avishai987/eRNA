# From https://umi-tools.readthedocs.io/en/latest/Single_cell_tutorial.html
# Using BWA for alignment
# Using featureCounts for assigning reads to genes

module snare_star:
    snakefile: "snare.smk"

use rule * from snare_star

hg_38_bwa_index = "/sci/data/reference_data/Homo_sapiens/Ensembl/GRCh38/Sequence/BWAIndex/genome.fa"

rule extract:
    input:
        script = "scripts/umi_tools/extract.sh",
        fastq1 = rules.snare_01_download_fastq.output.fastq1,
        fastq2 = rules.snare_01_download_fastq.output.fastq2,
        whitelist = rules.snare_02_download_counts.output.barcodes
    output:
        fastq1_extracted = "GSE126074_SNARE_seq/bwa/SRR8528318_1_extracted.fastq.gz",
        fastq2_extracted = "GSE126074_SNARE_seq/bwa/SRR8528318_2_extracted.fastq.gz"
    params:
        dir = "GSE126074_SNARE_seq/bwa"
    conda:
        "eRNA_umi_tools"
    log:
        "GSE126074_SNARE_seq/bwa/extract.log"
    shell:
        """
        mkdir -p {params.dir};
        {input.script} CCCCCCCCCCCCNNNNNNNNN {input.fastq1} {input.fastq2} {output.fastq1_extracted} \
        {output.fastq2_extracted} {input.whitelist} > {log} 2>&1
        """
rule align:
    input:
        script = "scripts/umi_tools/bwa.sh",
        fastq2 = rules.extract.output.fastq2_extracted
    output:
        bam = "GSE126074_SNARE_seq/bwa/SRR8528318.bam"
    params:
        dir = "GSE126074_SNARE_seq/bwa",
        cores = 4,
        ref_genome = hg_38_bwa_index
    conda:
        "eRNA_bwa"
    log:
        "GSE126074_SNARE_seq/bwa/bwa.log"
    shell:
        """
        mkdir -p {params.dir};
        {input.script} {params.cores} {params.ref_genome} {input.fastq2} {output.bam} > {log} 2>&1
        """

rule assign_to_genes:
    input:
        script = "scripts/umi_tools/featureCounts.sh",
        bam = rules.align.output.bam,
        gtf = rules.snare_09_bed_to_gtf_enhancers.output.gtf
    output:
        counts = "GSE126074_SNARE_seq/bwa/SRR8528318_counts.txt",
        sorted_tagged_bam = "GSE126074_SNARE_seq/bwa/SRR8528318_tagged_sorted.bam",
        bai_tagged_bam = "GSE126074_SNARE_seq/bwa/SRR8528318_tagged_sorted.bam.bai"
    params:
        dir = "GSE126074_SNARE_seq/bwa",
        feature_type = "regulatory_region",
        threads = 4
    conda:
        "eRNA_subread"
    log:
        "GSE126074_SNARE_seq/bwa/featureCounts.log"
    shell:
        """
        mkdir -p {params.dir};
        # remove XS tag from BAM (added by BWA), as it causes umi_tools to fail
        samtools view -h {input.bam} --remove-tag XS -b > {input.bam}.no_XS.bam;
        {input.script} {input.gtf} {input.bam}.no_XS.bam {output.counts} {params.threads} {params.feature_type} > {log} 2>&1;
        samtools sort {input.bam}.no_XS.bam.featureCounts.bam -o {output.sorted_tagged_bam} ;
        samtools index {output.sorted_tagged_bam};
        rm -f {input.bam}.no_XS.bam {input.bam}.featureCounts.bam;
        """

rule count_per_cell:
    input:
        sorted_tagged_bam = rules.assign_to_genes.output.sorted_tagged_bam
    output:
        count_matrix = "GSE126074_SNARE_seq/bwa/SRR8528318_counts_per_cell.txt",
    params:
        dir = "GSE126074_SNARE_seq/bwa"
    conda:
        "eRNA_umi_tools"
    log:
        "GSE126074_SNARE_seq/bwa/count_per_cell.log"
    shell:
        """
        umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XS --per-cell  -I {input.sorted_tagged_bam} -S {output.count_matrix} --wide-format-cell-counts > {log} 2>&1;
        """

