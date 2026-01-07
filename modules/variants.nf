#!/usr/bin/env nextflow

process CLAIR3 {

    //container 'docker.io/hkubal/clair3:latest'
    container 'docker.io/hkubal/clair3-gpu:latest'

    publishDir "${params.outdir}/03-variants", mode: 'copy'
    tag "${bam.simpleName}"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.variants.vcf"), path("${bam.simpleName}.variants.vcf.tbi")

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

    mv clair3_output/merge_output.vcf.gz ${bam.simpleName}.variants.vcf.gz
    mv clair3_output/merge_output.vcf.gz.tbi ${bam.simpleName}.variants.vcf.gz.tbi
    """ 
}

process VCF_STATS {
    
    container 'docker.io/staphb/bcftools:latest'
    publishDir "${params.outdir}/03-variants", mode: 'copy'
    tag "${vcf.simpleName}"

    input:
    tuple path(vcf), path(vcf_tbi)

    output:
    path("${vcf.simpleName}.query")

    script:
    """
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%FILTER\n' $vcf > ${vcf.simpleName}.query

    """
}
    
    