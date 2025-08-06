sbatch_params = "--mail-user=avishai.wizel@mail.huji.ac.il --export=ALL,conda_path=/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/bin/,lmod=/etc/profile.d/huji-lmod.sh"
log_dir = "/sci/labs/yotamd/lab_share/avishai.wizel/eRNA/PMID_33564857_breast_scRNA/logs/"

paths=" . /etc/profile.d/huji-lmod.sh; export conda_path=/sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/bin/ &&"
rule test:
    input:
        script = "./test/test.sh"
    output:
        output_name = "sss.txt"
    shell:
        '''
        {paths} {input.script} "my_log.log" "asd" "sdf" {output.output_name}
        '''

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
#                              GSE126074_SNARE_seq                             #
# ---------------------------------------------------------------------------- #

#scATAC-seq and scRNA-seq in the same cell

home_dir="GSE126074_SNARE_seq"
rule download_fastq:
    input:
        script = home_dir+"/scripts/download_fastq.sh"
    shell: 
        '''
        sbatch  {sbatch_params} \
        {input.script} \
        {home_dir}/01_raw_data/GSE126074/\
        {home_dir}\
        SRR8528318
        '''
#download available counts data
rule SNARE_seq_download_counts:
    shell: 
        '''
        cd {home_dir}/01_raw_data/GSE126074/
        wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074%5FCellLineMixture%5FSNAREseq%5FcDNA%5Fcounts.tsv.gz
        wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE126nnn/GSE126074/suppl/GSE126074%5FCellLineMixture%5FSNAREseq%5Fchromatin%5Fcounts.tsv.gz
        gzip -d GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv.gz
        gzip -d GSE126074_CellLineMixture_SNAREseq_chromatin_counts.tsv.gz
        head -n 1 GSE126074_CellLineMixture_SNAREseq_cDNA_counts.tsv | tr '\t' '\n' > barcodes.tsv
        '''

#align
rule star_alignment_SNARE_seq:
    input:
        script = home_dir +"/scripts/alignment.sh",
        fastq1 = home_dir + "/01_raw_data/GSE126074/SRR8528318_1.fastq",
        fastq2 = home_dir + "/01_raw_data/GSE126074/SRR8528318_2.fastq"

    shell: 
        '''
        sbatch  {sbatch_params} --mail-type=END\
        {input.script} \
        /sci/labs/yotamd/lab_share/avishai.wizel/python_envs/miniconda/envs/celescope/bin/STAR-avx2\
        PMID_33564857_breast_scRNA/01_CeleScope/hs_ensembl_99\
        {input.fastq1}\
        {input.fastq2}\
        {home_dir}/02_alignment/GSE126074_SNARE_seq_
        '''
# subset bam only for cells passed QC from paper
rule subset_BAM_SNARE_seq:
    input:
        script = home_dir +"/scripts/subset_bam.sh",
        bam_file = home_dir + "/02_alignment/GSE126074_SNARE_seq_Aligned.sortedByCoord.out.bam",
        filtered_barcode_list = home_dir + "/01_raw_data/GSE126074/barcodes.tsv",
    shell:
        '''
        sbatch  {sbatch_params}   --output={home_dir}/logs/subset_bam_%j.log  \
		{input.script} \
		{input.bam_file} \
        {input.filtered_barcode_list}\
        {home_dir}/03_count_erna/GSE126074_SNARE_seq_Aligned_filtered.sortedByCoord.out.bam
        '''
# convert bed to GTF- req for htseq
rule bed_to_gtf_SNARE_seq:
    input:
        script = "PMID_33564857_breast_scRNA/scripts/bed_to_gtf.sh",
        bed_file = home_dir+ "/03_count_erna/F03_nonExon_FANTOM5_typical_enhancers_no_prefix.bed"
    shell:
        '''
        sbatch  {sbatch_params} \
        {input.script}  {input.bed_file} {home_dir}/03_count_erna/F03_nonExon_FANTOM5_typical_enhancers_no_prefix.gtf
        '''
# for each eRNA locus, count reads with htseq
rule count_each_cell_SNARE_seq:
    input:
        script = "PMID_33564857_breast_scRNA/scripts/count_each_cell.sh",
        bam = home_dir + "/03_count_erna/GSE126074_SNARE_seq_Aligned_filtered.sortedByCoord.out.bam",
        gtf = home_dir + "/03_count_erna/F03_nonExon_FANTOM5_typical_enhancers_no_prefix.gtf"
    shell:
        '''
        sbatch  {sbatch_params} --output={home_dir}/logs/count_each_cell_%j.log \
		{input.script} \
		{input.bam} {input.gtf} {home_dir}/03_count_erna/GSE126074_SNARE_seq_FANTOM5_typical_enhancers_count_matrix.tsv
        '''
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
#find loci in atac-seq that intersect with enhancers annotation 
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