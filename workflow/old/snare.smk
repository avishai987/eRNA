



# ---------------------------------------------------------------------------- #
#                              GSE126074_SNARE_seq                             #
# ---------------------------------------------------------------------------- #

#snATAC-seq and snRNA-seq in the same cell


#-------------------------------- download data --------------------------------#

# 01- download raw fastq data
rule snare_01_download_fastq:
    input:
        script = "scripts/download_fastq.sh"
    output:
        log = "GSE126074_SNARE_seq/01_raw_data/GSE126074/log.txt",
        fastq1 = "GSE126074_SNARE_seq/01_raw_data/GSE126074/SRR8528318_1.fastq.gz",
        fastq2 = "GSE126074_SNARE_seq/01_raw_data/GSE126074/SRR8528318_2.fastq.gz"
    conda:
        "eRNA_sra-tools"
    shell: 
        '''
        mkdir -p GSE126074_SNARE_seq/01_raw_data/;
        {input.script} \
        GSE126074_SNARE_seq/01_raw_data/GSE126074/ \
        GSE126074_SNARE_seq \
        SRR8528318 \
        > {output.log} 2>&1;
        '''


#02- download available counts data
rule snare_02_download_counts:
    output:
        atac_counts = "GSE126074_SNARE_seq/01_raw_data/GSE126074/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv",
        rna_counts = "GSE126074_SNARE_seq/01_raw_data/GSE126074/GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv",
        barcodes = "GSE126074_SNARE_seq/01_raw_data/GSE126074/barcodes.tsv"
    
    shell: 
        '''
        cd GSE126074_SNARE_seq/01_raw_data/GSE126074/
        wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074%5FCellLineMixture%5FSNAREseq%5FcDNA%5Fcounts.tsv.gz
        wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074%5FCellLineMixture%5FSNAREseq%5Fchromatin%5Fcounts.tsv.gz
        gzip -d GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv.gz
        gzip -d GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv.gz
        head -n 1 GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv | tr '\t' '\n' > barcodes.tsv
        '''

#-------------------------------- process RNA-seq data --------------------------------#
#03- align fastq to genome with STAR
rule snare_03_alignment:
    input:
        script = "scripts/alignment.sh",
        fastq1 = rules.snare_01_download_fastq.output.fastq1,
        fastq2 = rules.snare_01_download_fastq.output.fastq2
    output:
        log = "GSE126074_SNARE_seq/02_alignment/log.txt",
        bam = "GSE126074_SNARE_seq/02_alignment/GSE126074_SNARE_seq_Aligned.sortedByCoord.out.bam"
    conda:
        "eRNA_star"
    shell: 
        '''
        mkdir -p GSE126074_SNARE_seq/02_alignment/; touch {output.log};
        {input.script} \
        {ensembl_113_star_index} \
        {input.fastq1} \
        {input.fastq2} \
        GSE126074_SNARE_seq/02_alignment/GSE126074_SNARE_seq_ \
        > GSE126074_SNARE_seq/02_alignment/log.txt 2>&1
        '''


#04- subset bam only for cells passed QC from paper 

rule snare_04_subset_BAM:
    input:
        script = "scripts/subset_bam.sh",
        bam_file =  rules.snare_03_alignment.output.bam,
        filtered_barcode_list =  rules.snare_02_download_counts.output.barcodes
    output:
        log = "GSE126074_SNARE_seq/03_subset_bam/log.txt",
        bam_filtered = "GSE126074_SNARE_seq/03_subset_bam/GSE126074_SNARE_seq_Aligned_filtered.sortedByCoord.out.bam"
    params:
        dir = "GSE126074_SNARE_seq/03_subset_bam/"
    conda:
        "eRNA_bamtools"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
		{input.bam_file} \
        {input.filtered_barcode_list}\
        {output.bam_filtered} \
        > {output.log} 2>&1
        '''
#-------------------------------- call peaks and count eRNA --------------------------------#

# First, we call peaks using MACS2 with the RNA reads from the filtered BAM file in order to identify transcribed regions.
# We than can intersect enhancers annotations with the peaks found in the RNA data to find eRNA  expression in out data.

# 05- call peaks with MACS2
rule snare_05_call_peaks:
    input:
        script = "scripts/masc2_run.sh",
        bam = rules.snare_04_subset_BAM.output.bam_filtered
    output:
        log = "GSE126074_SNARE_seq/04_peaks/log.txt",
        peaks = "GSE126074_SNARE_seq/04_peaks/SNARE_seq_peaks.narrowPeak"
    params:
        dir = "GSE126074_SNARE_seq/04_peaks/"
    conda:
        "eRNA_macs2"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
        {input.bam} \
        {params.dir} \
        "SNARE_seq" \
        > {output.log} 2>&1
        '''
# 06 - snare_06_get_genes_annotations
# get genes and regulatory features annotations from ensembl
# combine both annotations into one gff file (all_with_regulation_annotaions)
# all_with_regulation_annotaions contains all enhancers 
rule snare_06_get_genes_annotations:
    input:
        script = "scripts/get_genes_annotation.sh"
    output:
        log = "GSE126074_SNARE_seq/05_non_gene_body_peaks/download_log.txt",
        features_gff = "GSE126074_SNARE_seq/05_non_gene_body_peaks/Homo_sapiens.GRCh38.113.gff3",
        regulation_gff = "GSE126074_SNARE_seq/05_non_gene_body_peaks/Homo_sapiens.GRCh38.regulatory_features.v113.gff3",
        all_with_regulation_annotaions = "GSE126074_SNARE_seq/05_non_gene_body_peaks/Homo_sapiens.GRCh38.with_regulatory_features.v113.gff3",
    params:
        dir = "GSE126074_SNARE_seq/05_non_gene_body_peaks/"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
        {params.dir} \
        {output.all_with_regulation_annotaions} \
        > {output.log} 2>&1
        '''

# 07- get enhancers regions that not overlap with RNA-producing or gene-associated regions
rule snare_07_unique_enhancers:
    input:
        script = "scripts/filter_unique_enhancers_regions.sh",
        features_gff = rules.snare_06_get_genes_annotations.output.features_gff,
        regulation_gff = rules.snare_06_get_genes_annotations.output.regulation_gff
    output:
        log = "GSE126074_SNARE_seq/05_non_gene_body_peaks/unique_enhancers_log.txt",
        uniqe_enhancers = "GSE126074_SNARE_seq/05_non_gene_body_peaks/uniqe_enhancers.bed"
    conda: 
        "eRNA_bedtools"
    shell:
        '''
		{input.script} \
		{input.features_gff} \
		{input.regulation_gff} \
        {output.uniqe_enhancers} \
        > {output.log} 2>&1
        '''

# 08- count peaks overlapping with enhancers
# count peaks (narrowPeak) that overlap with :
# 1. unique enhancers regions (in file stats_file_enhancers)  
# 2.  with all genomic regions including non-unique enhancers (in file stats_file_all_regions)

rule snare_08_peaks_stats:
    input:
        script = "scripts/count_peaks_for_all_annotations.sh",
        peaks = rules.snare_05_call_peaks.output.peaks,
        gff = rules.snare_06_get_genes_annotations.output.all_with_regulation_annotaions,
        uniqe_enhancers = rules.snare_07_unique_enhancers.output.uniqe_enhancers
    output:
        log = "GSE126074_SNARE_seq/05_non_gene_body_peaks/stats_log.txt",
        stats_file_all_regions = "GSE126074_SNARE_seq/05_non_gene_body_peaks/stats_all_regions.tsv",
        stats_file_enhancers= "GSE126074_SNARE_seq/05_non_gene_body_peaks/stats_unique_enhancers.tsv"
    conda: 
        "eRNA_bedtools"
    params:
        dir = "GSE126074_SNARE_seq/05_non_gene_body_peaks/"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
		{input.peaks} \
		{input.gff} \
        14 \
        {output.stats_file_all_regions} \
        > {output.log} 2>&1; 
	    {input.script} \
		{input.peaks} \
		{input.uniqe_enhancers} \
        15 \
        {output.stats_file_enhancers} \
        > {output.log} 2>&1;
        '''

#-------------------------------- count eRNA --------------------------------#
# convert bed to GTF- req for htseq
rule snare_09_bed_to_gtf_enhancers:
    input:
        script = "scripts/bed_to_gtf.sh",
        bed_file =rules.snare_07_unique_enhancers.output.uniqe_enhancers,
    output:
        log = "GSE126074_SNARE_seq/06_count_eRNA/bed_to_gtf.log",
        gtf = "GSE126074_SNARE_seq/06_count_eRNA/regulatory_ENCODE_cCREs_GRCh38_.gtf"
    params:
        dir = "GSE126074_SNARE_seq/06_count_eRNA/"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
        {input.bed_file} {output.gtf}\
        > {output.log} 2>&1;
        '''
# for each eRNA locus, count reads with htseq
rule snare_10_count_each_cell:
    input:
        script = "scripts/count_each_cell.sh",
        bam = rules.snare_04_subset_BAM.output.bam_filtered,
        gtf = rules.snare_09_bed_to_gtf_enhancers.output.gtf
    output:
        log = "GSE126074_SNARE_seq/06_count_eRNA/count_each_cell.log",
        mat = "GSE126074_SNARE_seq/06_count_eRNA/erna_exrs_matrix.tsv"
    params:
        dir = "GSE126074_SNARE_seq/06_count_eRNA/"
    conda:
        "eRNA_htseq"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
		{input.bam} {input.gtf} {output.mat}\
        > {output.log} 2>&1;
        '''
#-------------------------------- associate eRNA with chromatin regions --------------------------------#
# convert chromatin counts regions (first column) to bed format
rule snare_11_chromatin_regions_to_bed:
    input:
        script = "scripts/chromatin_regions_to_bed.sh",
        chromatin_file = rules.snare_02_download_counts.output.atac_counts
    output:
        log = "GSE126074_SNARE_seq/07_regions_to_erna/chromatin_regions_to_bed.log",
        bed_file = "GSE126074_SNARE_seq/07_regions_to_erna/chromatin_regions.bed"
    params:
        dir = "GSE126074_SNARE_seq/07_regions_to_erna/"
    shell:
        '''
        mkdir -p {params.dir};
        {input.script} \
        {input.chromatin_file} {output.bed_file}\
        > {output.log} 2>&1;
        '''

# associate enhancers with chromatin
rule snare_12_associate_enhancers_with_regions:
    input:
        script = "scripts/associate_inhancers_with_regions.sh",
        regions_file = rules.snare_11_chromatin_regions_to_bed.output.bed_file,
        enhancers_file = rules.snare_07_unique_enhancers.output.uniqe_enhancers
    output:
        log = "GSE126074_SNARE_seq/07_regions_to_erna/associate_enhancers_with_regions.log",
        associated_file = "GSE126074_SNARE_seq/07_regions_to_erna/associated_enhancers_with_chromatin.bed"
    shell:
        '''
        {input.script} \
        {input.regions_file} \
        {input.enhancers_file} \
        {output.associated_file} \
        > {output.log} 2>&1;
        '''
# analyze counts with GSE126074_SNARE_seq/08_counts_analysis/analysis.ipynb
rule snare_13_counts_analysis:
        input:
            script = "GSE126074_SNARE_seq/08_counts_analysis/analysis.ipynb",
            counts = rules.snare_10_count_each_cell.output.mat,
            associated_file = rules.snare_12_associate_enhancers_with_regions.output.associated_file
        output:
            report =  'GSE126074_SNARE_seq/08_counts_analysis/report/analysis.html'
        params:
            dir = "GSE126074_SNARE_seq/08_counts_analysis/report/"
        conda:
            "eRNA_jupyter"
        shell:
            '''

            mkdir -p GSE126074_SNARE_seq/08_counts_analysis/report;
            jupyter nbconvert --to html --execute {input.script} --output-dir {params.dir} --output analysis.html  --log-level=INFO
            '''
    




       
rule snare_all:
    input:
        rules.snare_01_download_fastq.output.log,
        # rules.snare_02_download_counts.output.log,
        rules.snare_03_alignment.output.log,
        rules.snare_04_subset_BAM.output.log,
        rules.snare_05_call_peaks.output.log,
        rules.snare_06_get_genes_annotations.output.log,
        rules.snare_07_unique_enhancers.output.log,
        rules.snare_08_peaks_stats.output.log,
        rules.snare_09_bed_to_gtf_enhancers.output.log,
        rules.snare_10_count_each_cell.output.log,
        rules.snare_11_chromatin_regions_to_bed.output.log,
        rules.snare_12_associate_enhancers_with_regions.output.log,
        rules.snare_13_counts_analysis.output.report
    output:
        "./GSE126074_SNARE_seq/snare_all_done.txt"
    shell:
        "echo 'done' > {output}"



