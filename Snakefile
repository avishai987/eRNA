sbatch_params = "--mail-user=avishai.wizel@mail.huji.ac.il --export=ALL,conda_path=/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/bin/,lmod=/etc/profile.d/huji-lmod.sh"
log_dir = "/sci/labs/yotamd/lab_share/avishai.wizel/eRNA/PMID_33564857_breast_scRNA/logs/"

paths=" . /etc/profile.d/huji-lmod.sh; export conda_path=/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/bin/ &&"


rule test:
    input:
        script = "./test/test.sh"
    output:
        log = "test_log.txt"

    conda:
        "eRNA"
    shell:
        '''
        {input.script}  > {output.log} 2>&1
        '''



# ---------------------------------------------------------------------------- #
#                              GSE126074_SNARE_seq                             #
# ---------------------------------------------------------------------------- #

#snATAC-seq and snRNA-seq in the same cell
#-------------------------------- download data --------------------------------#
home_dir="GSE126074_SNARE_seq"
rule snare_download_fastq:
    input:
        script = "scripts/download_fastq.sh"
    output:
        log = "GSE126074_SNARE_seq/01_raw_data/GSE126074/log.txt",
        fastq1 = "GSE126074_SNARE_seq/01_raw_data/GSE126074/SRR8528318_1.fastq.gz",
        fastq2 = "GSE126074_SNARE_seq/01_raw_data/GSE126074/SRR8528318_2.fastq.gz"
    conda:
        "eRNA"
    shell: 
        '''
        mkdir -p GSE126074_SNARE_seq/01_raw_data/;
        {input.script} \
        GSE126074_SNARE_seq/01_raw_data/GSE126074/ \
        GSE126074_SNARE_seq \
        SRR8528318 \
        > {output.log} 2>&1;
        '''


#download available counts data
rule SNARE_seq_download_counts:
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
#align fastq to genome with STAR
rule snare_star_alignment:
    input:
        script = "scripts/alignment.sh",
        fastq1 = rules.snare_download_fastq.output.fastq1,
        fastq2 = rules.snare_download_fastq.output.fastq2
    output:
        log = "GSE126074_SNARE_seq/02_alignment/log.txt",
        bam = "GSE126074_SNARE_seq/02_alignment/GSE126074_SNARE_seq_Aligned.sortedByCoord.out.bam"
    conda:
        "eRNA"
    shell: 
        '''
        mkdir -p GSE126074_SNARE_seq/02_alignment/; touch {output.log};
        {input.script} \
        /sci/data/reference_data/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex/Ens113/ \
        {input.fastq1} \
        {input.fastq2} \
        GSE126074_SNARE_seq/02_alignment/GSE126074_SNARE_seq_ \
        > GSE126074_SNARE_seq/02_alignment/log.txt 2>&1
        '''


# subset bam only for cells passed QC from paper 

rule subset_BAM_SNARE_seq:
    input:
        script = "scripts/subset_bam.sh",
        bam_file =  rules.snare_star_alignment.output.bam,
        filtered_barcode_list =  rules.SNARE_seq_download_counts.output.barcodes
    output:
        log = "GSE126074_SNARE_seq/03_subset_bam/log.txt",
        bam_filtered = "GSE126074_SNARE_seq/03_subset_bam/GSE126074_SNARE_seq_Aligned_filtered.sortedByCoord.out.bam"
    params:
        dir = "GSE126074_SNARE_seq/03_subset_bam/"
    conda:
        "eRNA"
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

rule call_peaks_snare:
    input:
        script = "scripts/masc2_run.sh",
        bam = rules.subset_BAM_SNARE_seq.output.bam_filtered
    output:
        log = "GSE126074_SNARE_seq/04_peaks/log.txt",
        peaks = "GSE126074_SNARE_seq/04_peaks/SNARE_seq_peaks.narrowPeak"
    params:
        dir = "GSE126074_SNARE_seq/04_peaks/"
    conda:
        "eRNA"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
        {input.bam} \
        {params.dir} \
        "SNARE_seq" \
        > {output.log} 2>&1
        '''
# get genes and regulatory features annotations from ensembl
# combine both annotations into one gff file (all_with_regulation_annotaions)
# all_with_regulation_annotaions contains all enhancers 
rule get_genes_annotations:
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

# get enhancers regions that not overlap with RNA-producing or gene-associated regions 
rule get_unique_enhancers_regions:
    input:
        script = "scripts/filter_unique_enhancers_regions.sh",
        features_gff = rules.get_genes_annotations.output.features_gff,
        regulation_gff = rules.get_genes_annotations.output.regulation_gff,
    output:
        log = protected("GSE126074_SNARE_seq/05_non_gene_body_peaks/unique_enhancers_log.txt"),
        uniqe_enhancers = "GSE126074_SNARE_seq/05_non_gene_body_peaks/uniqe_enhancers.bed"
    shell:
        '''
		{input.script} \
		{input.features_gff} \
		{input.regulation_gff} \
        {output.uniqe_enhancers} \
        > {output.log} 2>&1
        '''

# count peaks (narrowPeak) that overlap with :
# 1. unique enhancers regions (in file stats_file_enhancers)  
# 2.  with all genomic regions including non-unique enhancers (in file stats_file_all_regions)

rule snare_peaks_stats:
    input:
        script = "scripts/count_peaks_for_all_annotations.sh",
        peaks = rules.call_peaks_snare.output.peaks,
        gff = rules.get_genes_annotations.output.all_with_regulation_annotaions,
        uniqe_enhancers = rules.get_unique_enhancers_regions.output.uniqe_enhancers
    output:
        log = protected("GSE126074_SNARE_seq/05_non_gene_body_peaks/stats_log.txt"),
        stats_file_all_regions = "GSE126074_SNARE_seq/05_non_gene_body_peaks/stats_all_regions.tsv",
        stats_file_enhancers= "GSE126074_SNARE_seq/05_non_gene_body_peaks/stats_unique_enhancers.tsv"
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
rule snare_bed_to_gtf_enhancers:
    input:
        script = "scripts/bed_to_gtf.sh",
        bed_file =rules.get_unique_enhancers_regions.output.uniqe_enhancers,
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
rule snare_count_each_cell:
    input:
        script = "scripts/count_each_cell.sh",
        bam = rules.subset_BAM_SNARE_seq.output.bam_filtered,
        gtf = rules.snare_bed_to_gtf_enhancers.output.gtf
    output:
        log = "GSE126074_SNARE_seq/06_count_eRNA/count_each_cell.log",
        mat = "GSE126074_SNARE_seq/06_count_eRNA/erna_exrs_matrix.tsv"
    shell:
        '''
		{input.script} \
		{input.bam} {input.gtf} {output.mat}\
        > {output.log} 2>&1;
        '''
#-------------------------------- associate eRNA with chromatin regions --------------------------------#
# convert chromatin counts regions (first column) to bed format
rule snare_chromatin_regions_to_bed:
    input:
        script = "scripts/chromatin_regions_to_bed.sh",
        chromatin_file = rules.SNARE_seq_download_counts.output.atac_counts
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
rule snare_associate_enhancers_with_regions:
    input:
        script = "scripts/associate_inhancers_with_regions.sh",
        regions_file = rules.snare_chromatin_regions_to_bed.output.bed_file,
        enhancers_file = rules.get_unique_enhancers_regions.output.uniqe_enhancers
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
#-------------------------------- with known specific cell type enhancers --------------------------------#

# download enhancers annotations by cell type
rule download_annotations:
    shell:
        '''
        sbatch  {sbatch_params}   --output={home_dir}/logs/{rule}_%j.log  \
	      --wrap=' cd {home_dir}/01_raw_data/enhancers_anotation/;
        wget http://www.enhanceratlas.org/data/download/enhancer/hs/H1.bed;
        wget http://www.enhanceratlas.org/data/download/enhancer/hs/K562.bed;
        wget http://www.enhanceratlas.org/data/download/enhancer/hs/GM12878.bed;
        wget http://www.enhanceratlas.org/data/download/enhancer/hs/BJ.bed'
        '''
# intersect enhancers annotation
rule intersect_annotations:
    input:
        script = home_dir +"/scripts/intersect_annotations.sh",
        bed_1 = home_dir + "/01_raw_data/enhancers_anotation/BJ.bed",
        bed_2 = home_dir + "/01_raw_data/enhancers_anotation/GM12878.bed",
        bed_3 = home_dir + "/01_raw_data/enhancers_anotation/H1.bed",
        bed_4 = home_dir + "/01_raw_data/enhancers_anotation/K562.bed"
    params:
        out_file = home_dir + "/04_annotations_intersect/intersected_annotation.bed"  
    shell:
        '''
        sbatch  {sbatch_params}   --output={home_dir}/logs/{rule}_%j.log  \
	    {input.script} {input.bed_1} {input.bed_2} {input.bed_3} {input.bed_4} {params.out_file}
       '''
#find regions in atac-seq that intersect with enhancers annotation 
rule filtered_bed_for_atac:
    input:
        script = home_dir +"/scripts/filter_atac_loci.sh",
        atac_seq_file = home_dir + "/01_raw_data/GSE126074/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv",
        annotations_bed_file = rules.intersect_annotations.params.out_file
    params:
        out_file = home_dir + "/04_annotations_intersect/filtered_atac_seq.tsv"  
    shell:
        '''
        sbatch  {sbatch_params}   --output={home_dir}/logs/{rule}_%j.log  \
	    {input.script} {input.atac_seq_file} {input.annotations_bed_file} {params.out_file}
       '''
# subset sc-ATACseq to enhancers loci and scRNA-seq to var genes
rule filter_rna_and_atac:
      input:
        script = home_dir +"/scripts/run_Rscript.sh",
        r_scrpit = home_dir +"/r_notebooks/filter_data.Rmd"
      params:
        report =  home_dir + '/05_filtered_data/filter_data.html'
      shell:
        '''
        sbatch {sbatch_params} --mem=8G --time=1:0:0 --output={home_dir}/logs/{rule}_%j.log  \
	      {input.script} \
	      {input.r_scrpit} {params.report}
       '''
       
rule snare_all:
    input:
        rules.snare_download_fastq.output.log,
        # rules.SNARE_seq_download_counts.output.log,
        rules.snare_star_alignment.output.log,
        rules.subset_BAM_SNARE_seq.output.log,
        rules.call_peaks_snare.output.log,
        rules.get_genes_annotations.output.log,
        rules.get_unique_enhancers_regions.output.log,
        rules.snare_peaks_stats.output.log,
        rules.snare_bed_to_gtf_enhancers.output.log,
        rules.snare_count_each_cell.output.log
    shell:
        "echo 'done'"
# ---------------------------------------------------------------------------- #
#                          PMID_33564857_breast_scRNA                          #
# ---------------------------------------------------------------------------- #

rule index_bam:
    input:
        script = "PMID_33564857_breast_scRNA/scripts/index_bam.sh",
        bam_file = "/sci/labs/yotamd/lab_share/avishai.wizel/eRNA/PMID_33564857_breast_scRNA/01_CeleScope/breast_sample/outs/breast_sample_Aligned.sortedByCoord.out.bam"
    shell:
        "sbatch --chdir /sci/labs/yotamd/lab_share/avishai.wizel/eRNA/PMID_33564857_breast_scRNA/01_CeleScope/breast_sample/outs/ {sbatch_params} \
 {input.script} {input.bam_file} " 
 
 # subset the BAM to include reads from cells that passed star filteration
rule subset_BAM:
    input:
        script = "PMID_33564857_breast_scRNA/scripts/subset_BAM.sh",
        bam_file = "PMID_33564857_breast_scRNA/01_CeleScope/breast_sample/outs/breast_sample_Aligned.sortedByCoord.out.bam",
        filtered_barcode_list = "PMID_33564857_breast_scRNA/01_CeleScope/breast_sample/outs/filtered/barcodes.tsv"
    shell:
        "sbatch  {sbatch_params}   --output='PMID_33564857_breast_scRNA/logs/subset_bam_%j.log' --mail-type=NONE  \
		{input.script} \
		'PMID_33564857_breast_scRNA/05_erna_count/' \
		{input.bam_file} \
        {input.filtered_barcode_list}" 
         
# index the subset of BAM
rule index_subset_bam:
    input:
        script = "PMID_33564857_breast_scRNA/scripts/index_bam_2.sh",
        bam_file = 'PMID_33564857_breast_scRNA/05_erna_count/filtered_bam.bam' 
    shell:
        "sbatch  {sbatch_params} \
	 {input.script}  {input.bam_file} "

# convert bed to GTF- req for htseq
rule bed_to_gtf:
    input:
        script = "PMID_33564857_breast_scRNA/scripts/bed_to_gtf.sh",
        bed_file = "PMID_33564857_breast_scRNA/04_eRNA_loci_intersection/PMID_33564857_erna_peaks_non_exons.bed"
    shell:
        "sbatch  {sbatch_params} \
	 {input.script}  {input.bed_file} 'PMID_33564857_breast_scRNA/05_erna_count/PMID_33564857_erna_peaks_non_exons.gtf' "

# for each eRNA locus, count reads with htseq
rule count_each_cell:
    input:
        script = "PMID_33564857_breast_scRNA/scripts/count_each_cell.sh",
        bam = "PMID_33564857_breast_scRNA/01_CeleScope/breast_sample/outs/breast_sample_Aligned.sortedByCoord.out.bam",
        gtf = "PMID_33564857_breast_scRNA/05_erna_count/PMID_33564857_erna_peaks_non_exons.gtf"
    shell:
        "sbatch  {sbatch_params} \
		{input.script} \
		{input.bam} {input.gtf} PMID_33564857_erna_peaks_non_exons_count_matrix.txt" 



# ---------------------------------------------------------------------------- #
#                              GSE126074_SNARE_seq-  BABEL                     #
# ---------------------------------------------------------------------------- #
# use BABEL to predict chromatin accessability 

# create input
rule create_h5_files:
    input:
        script = home_dir +"/scripts/run_tsv_to_h5.sh",
        py_script = home_dir + "/scripts/tsv_to_h5.py",
        path_to_rna_tsv = home_dir + "/01_raw_data/GSE126074/GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv",
        path_to_atac_tsv =home_dir + "/01_raw_data/GSE126074/GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv"
    shell:
        '''
        sbatch  {sbatch_params}   --output={home_dir}/logs/create_h5_files_%j.log  \
	    {input.script} {input.py_script} {input.path_to_rna_tsv} {input.path_to_atac_tsv}  {home_dir}/BABEL/01_h5_files/GSE126074_SNARE_seq.h5
        '''
        
        

       
rule run_babel_SNARE_seq:
    input:
        script = home_dir +"/scripts/run_babel.sh",
        py_script = home_dir + "/BABEL/babel_repo/bin/train_model.py",
        path_to_rna_atac_h5 = home_dir + "/BABEL/01_h5_files/GSE126074_SNARE_seq.h5"
    params:
        out_dir = home_dir + "/BABEL/02_run/"
    shell:
        '''
        sbatch  {sbatch_params}  --output={home_dir}/logs/{rule}_%j.log  \
	    {input.script} {input.py_script} {input.path_to_rna_atac_h5}   {input.params.out_dir}
        '''


# ---------------------------------------------------------------------------- #
#                               10X_PBMC_GEX_ATAC                              #
# ---------------------------------------------------------------------------- #

rule download_10X_atac:
    #data description link: https://support.10xgenomics.com/single-cell-multiome-atac-gex/datasets/1.0.0/pbmc_granulocyte_sorted_10k
    shell:
        '''
        mkdir -p ./10X_PBMC/01_raw_data
        cd ./10X_PBMC/01_raw_data
        wget https://github.com/wukevin/babel/raw/refs/heads/main/data/10x/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5
        '''