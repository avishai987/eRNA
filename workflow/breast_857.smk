
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

