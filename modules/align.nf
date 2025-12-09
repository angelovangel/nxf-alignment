#!/usr/bin/env nextflow

process DORADO_ALIGN {

    container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/01-align", mode: 'copy', pattern: '*{bam,bai}'
    tag "${reads.simpleName}"

    input:
        tuple path(ref), path(reads)

    output:
        path("*.{bam,bai}")

    script:
    """
    dorado aligner ${ref} ${reads} | samtools sort -o ${reads.simpleName}.align.bam
    samtools index ${reads.simpleName}.align.bam
    """
}

process MAKE_BEDFILE {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    
    input:
        path ref

    output:
        path "fallback.bed"

    script:
    """
    samtools faidx ${ref}
    awk -v OFS='\t' '{print \$1, 0, \$2, \$1}' ${ref}.fai > fallback.bed 
    """
}

process BEDTOOLS_COV {
    container 'docker.io/biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*hist.tsv'
    tag "${bam.simpleName}"

    input:
        tuple path(bam), path(bai), path(bed)

    output:
        path "*hist.tsv"

    script:
    """
    echo -e "chr\tstart\tend\tlabel\tdepth\tbases_at_depth\tsize\tpercent_at_depth" > ${bam.simpleName}.hist.tsv
    bedtools coverage -a ${bed} -b ${bam} -hist >> ${bam.simpleName}.hist.tsv
    """
}