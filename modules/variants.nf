#!/usr/bin/env nextflow

process VCF_CLAIR3 {

    //container 'docker.io/hkubal/clair3:latest'
    container 'docker.io/hkubal/clair3-gpu:latest'

    publishDir "${params.outdir}/03-variants", mode: 'copy'
    errorStrategy 'ignore'
    tag "${bam.simpleName}"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.vcf.gz"), path("${bam.simpleName}.vcf.gz.tbi")

    script:
    def model = "${params.clair3_model}"
    def platform = "${params.clair3_platform}"
    """
    samtools faidx $ref

    /opt/bin/run_clair3.sh \
    --bam_fn=$bam \
    --ref_fn=$ref \
    --bed_fn=$bedfile \
    --platform=$platform \
    --model_path="/opt/models/$model" \
    --threads=8 \
    --output="clair3_output"

    mv clair3_output/merge_output.vcf.gz ${bam.simpleName}.vcf.gz
    mv clair3_output/merge_output.vcf.gz.tbi ${bam.simpleName}.vcf.gz.tbi
    """ 
}

process VCF_DEEPVARIANT {
    container 'docker.io/google/deepvariant:1.10.0-beta-gpu'
    publishDir "${params.outdir}/03-variants", mode: 'copy'
    tag "${bam.simpleName}"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.vcf.gz"), path("${bam.simpleName}.vcf.gz.tbi")

    script:
    """
    samtools faidx $ref

    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${params.deepvariant_model} \
    --ref=$ref \
    --reads=$bam \
    --output_vcf=${bam.simpleName}.vcf.gz
    --num_shards=8
    """ 
}

process VCF_STATS {
    
    container 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'
    //publishDir "${params.outdir}/03-variants", mode: 'copy'
    errorStrategy 'ignore'
    tag "${vcf.simpleName}"

    input:
    tuple path(vcf), path(vcf_tbi)

    output:
    path("${vcf.simpleName}.query")

    script:
    """
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%FILTER\t%QUAL\n' $vcf > ${vcf.simpleName}.query

    """
}

process VCF_ANNOTATE {
    container 'docker.io/dceoy/snpeff:latest'
    containerOptions '--entrypoint ""' 
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/04-annotations", mode: 'copy'
    tag "${vcf.simpleName}"

    input:
    tuple path(vcf), path(vcf_tbi)

    output:
    path("${vcf.simpleName}.ann.vcf")
    path("${vcf.simpleName}.stats.csv"), emit: ch_vcfann_stats

    script:
    """
    # download snpEff database here
    mkdir -p ./snpeff_data
    snpEff ann -dataDir \$PWD/snpeff_data -csvStats ${vcf.simpleName}.stats.csv ${params.anno_db} $vcf > ${vcf.simpleName}.ann.vcf
    
    """
}

process VCF_ANNOTATE_REPORT {
    publishDir "${params.outdir}/04-annotations", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path(ann_stats)

    output:
    path("variants_annotation_report.html")

    script:
    """
    make-variants-report.py $ann_stats -o variants_annotation_report.html
    """
}
    