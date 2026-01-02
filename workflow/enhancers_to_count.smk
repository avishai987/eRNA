#  get_genes_annotations 
# get genes and regulatory features annotations from ensembl
# combine both annotations into one gff file (all_with_regulation_annotaions)
# all_with_regulation_annotaions contains all enhancers 
rule get_genes_annotations:
    input:
        script = "scripts/get_genes_annotation.sh"
    output:
        log = "Analysis/enhancers/ensembl/raw/download_log.txt",
        features_gff = "Analysis/enhancers/ensembl/raw/Homo_sapiens.GRCh38.113.gff3",
        regulation_gff = "Analysis/enhancers/ensembl/raw/Homo_sapiens.GRCh38.regulatory_features.v113.gff3",
        all_with_regulation_annotaions = "Analysis/enhancers/ensembl/raw/Homo_sapiens.GRCh38.with_regulatory_features.v113.gff3",
    params:
        dir = "Analysis/enhancers/ensembl/raw/"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
        {params.dir} \
        {output.all_with_regulation_annotaions} \
        > {output.log} 2>&1
        '''

# get enhancers regions that not overlap with RNA-producing or gene-associated regions
rule get_unique_enhancers:
    input:
        script = "scripts/filter_unique_enhancers_regions.sh",
        features_gff = rules.get_genes_annotations.output.features_gff,
        regulation_gff = rules.get_genes_annotations.output.regulation_gff
    output:
        log = "Analysis/enhancers/ensembl/unique_enhancers_log.txt",
        uniqe_enhancers = config["enhancers_to_count"]["encode_bed"]
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

# convert bed to GTF- req for htseq
rule bed_to_gtf_enhancers:
    input:
        script = "scripts/bed_to_gtf.sh",
        bed_file =rules.get_unique_enhancers.output.uniqe_enhancers,
    output:
        log = "Analysis/enhancers/ensembl/bed_to_gtf.log",
        gtf = config["enhancers_to_count"]["encode_gtf"]
    params:
        dir = "Analysis/enhancers/ensembl/"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
        {input.bed_file} {output.gtf}\
        > {output.log} 2>&1;
        '''
rule get_enhancers_meta:
    input:
        script = "scripts/get_enhancers_metadata.sh",
        uniqe_enhancers = rules.get_unique_enhancers.output.uniqe_enhancers,
        features_gff = rules.get_genes_annotations.output.features_gff
    output:
        log = "Analysis/enhancers/ensembl/enhancers_metadata_log.txt",
        metadata = "Analysis/enhancers/ensembl/ensembl_enhancers_metadata.txt"
    params:
        dir = "Analysis/enhancers/ensembl/"
    shell:
        '''
        mkdir -p {params.dir};
		{input.script} \
        {input.uniqe_enhancers} \
        {input.features_gff} \
        {output.metadata} \
        > {output.log} 2>&1;
        '''