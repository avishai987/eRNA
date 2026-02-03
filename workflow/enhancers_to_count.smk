#  get_genes_annotations 
# get genes annotations from ensembl
# download genes annotation https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz
rule download_genes_annotation:
    output:
        genes_annotation = "Analysis/enhancers/genes_annotation/Homo_sapiens.GRCh38.113.gff3",
        hg19ToHg38_over_chain = "Analysis/enhancers/genes_annotation/hg19ToHg38.over.chain.gz"
    params:
        dir = "Analysis/enhancers/genes_annotation/"
    shell:
        '''
        wget https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz -P {params.dir}
        gunzip {params.dir}Homo_sapiens.GRCh38.113.gff3.gz
        wget https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz -P {params.dir}
        '''

#download regulatory features and convert to bed
rule set_ensemble_regulatory_features:
    input:
        script = "scripts/ensembl_gff_to_bed.sh",
    output:
        regulation_gff = "Analysis/enhancers/ensembl/Homo_sapiens.GRCh38.regulatory_features.v113.gff3",
        regulation_bed = "Analysis/enhancers/ensembl/Homo_sapiens.GRCh38.regulatory_features.v113.bed"
    params:
        dir = "Analysis/enhancers/ensembl/"
    shell:
        '''
        mkdir -p {params.dir};
        wget -O {output.regulation_gff}.gz https://ftp.ensembl.org/pub/release-113/regulation/homo_sapiens/GRCh38/annotation/Homo_sapiens.GRCh38.regulatory_features.v113.gff3.gz
        gunzip {output.regulation_gff}.gz
        {input.script} {output.regulation_gff} {output.regulation_bed}
        '''

#from https://bioinformatics.mdanderson.org/Supplements/Super_Enhancer/TCEA_website/parts/2_Samples_and_annotation.html
#download, ungzip, convert to hg38 and to GTF, remove "chr" from first col
rule set_tcea_fantom_enhancers:
    input:
        hg19ToHg38_over_chain = rules.download_genes_annotation.output.hg19ToHg38_over_chain
    output:
        enhancers_bed = "Analysis/enhancers/tcea_FANTOM/F03_nonExon_FANTOM5_typical_enhancers_hg38.bed"
    params:
        dir = "Analysis/enhancers/tcea_FANTOM/"
    conda: 
        "eRNA_bedtools"
    shell:
        '''
        wget https://bioinformatics.mdanderson.org/Supplements/Super_Enhancer/7_Typical_enhancer_Dec2020_update_FANTOM_60k/F03_nonExon_FANTOM5_typical_enhancers.bed.gz -P {params.dir}
        gunzip {params.dir}F03_nonExon_FANTOM5_typical_enhancers.bed.gz
        liftOver {params.dir}F03_nonExon_FANTOM5_typical_enhancers.bed {input.hg19ToHg38_over_chain} {output.enhancers_bed} {params.dir}unMapped.bed
        # remove "chr" from first column
        sed -i 's/chr//g' {output.enhancers_bed}
        # keep only standard chromosomes
        grep -E '^[0-9XYM]{{1,2}}[[:space:]]' {output.enhancers_bed} > {output.enhancers_bed}.tmp
        mv {output.enhancers_bed}.tmp {output.enhancers_bed}
        '''   

rule set_tcea_super_enhancers:
    input:
        hg19ToHg38_over_chain = rules.download_genes_annotation.output.hg19ToHg38_over_chain,
        script = "scripts/set_tcea_super_enhancers.sh"
    output:
        super_enhancers_bed = "Analysis/enhancers/tcea_super_enhancers/tcea_superEnhancers.bed"
    params:
        dir = "Analysis/enhancers/tcea_super_enhancers/"
    conda: 
        "eRNA_bedtools"
    shell:
        '''
        mkdir -p {params.dir};
        {input.script} {params.dir} {input.hg19ToHg38_over_chain} {output.super_enhancers_bed}
        '''

ENHANCER_BED_MAP = {
    "ensembl": rules.set_ensemble_regulatory_features.output.regulation_bed,
    "tcea_FANTOM": rules.set_tcea_fantom_enhancers.output.enhancers_bed,
    "tcea_super_enhancers": rules.set_tcea_super_enhancers.output.super_enhancers_bed
}

rule filter_all_bed:
    input:
        expand("Analysis/enhancers/{name}/{name}_filtered_enhancers.bed", name=ENHANCER_BED_MAP.keys())

rule filter_bed:
    input:
        enhancer = lambda wildcards: ENHANCER_BED_MAP[wildcards.name],
        genes_annotation = rules.download_genes_annotation.output.genes_annotation,
        bed_to_gtf_script = "scripts/bed_to_gtf.sh"
    output:
        log = "Analysis/enhancers/{name}/filtered_enhancers_log.txt",
        filtered_enhancers = "Analysis/enhancers/{name}/{name}_filtered_enhancers.bed",
        filtered_enhancers_gtf = "Analysis/enhancers/{name}/{name}_filtered_enhancers.gtf"
    conda: 
        "eRNA_bedtools"
    shell:
        """
        bedtools intersect -a {input.enhancer} -b <(awk '$3 == "exon"' {input.genes_annotation}) -v  > {output.filtered_enhancers}  2> {output.log}
        {input.bed_to_gtf_script} {output.filtered_enhancers} {output.filtered_enhancers_gtf} 2>&1 | tee -a {output.log}
        """

ENHANCER_metadata_MAP = {
    "ensembl": rules.set_ensemble_regulatory_features.output.regulation_bed,
    "tcea_FANTOM": rules.set_tcea_fantom_enhancers.output.enhancers_bed,
    "tcea_super_enhancers": rules.set_tcea_super_enhancers.output.super_enhancers_bed
}

rule get_all_enhancers_meta:
    input:
        expand("Analysis/enhancers/{name}/{name}_enhancers_metadata.txt", name=ENHANCER_metadata_MAP.keys())
        
rule get_enhancers_meta:
    input:
        script = "scripts/get_enhancers_metadata.sh",
        enhancers_bed = lambda wildcards: ENHANCER_metadata_MAP[wildcards.name],
        genes_annotation = rules.download_genes_annotation.output.genes_annotation
    output:
        log = "Analysis/enhancers/{name}/{name}_enhancers_metadata_log.txt",
        metadata = "Analysis/enhancers/{name}/{name}_enhancers_metadata.txt"
    params:
        dir = "Analysis/enhancers/{name}/"
    shell:
        '''
        mkdir -p {params.dir};
        {input.script} \
        {input.enhancers_bed} \
        {input.genes_annotation} \
        {output.metadata} \
        > {output.log} 2>&1;
        '''
