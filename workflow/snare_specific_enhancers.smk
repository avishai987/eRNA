# workflow for measurement of snare-seq specific cell type enhancers
module snare_all_enhancers:
    snakefile: "Snakefile"

use rule * from snare_all_enhancers 

rule test_specific_enhancers:
    input:
        rules.snare_download_fastq.output.fastq1,
        rules.snare_download_fastq.output.fastq2
    output:
        "test_output.txt"
    shell:
        """
        echo "Test successful!" > {output}
        """

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