#!/usr/bin/env nextflow

process MERGE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    errorStrategy 'ignore' //because some barcodes defined in the samplesheet might be missing in the data
    tag "${barcode} == ${samplename}"

    publishDir "$params.outdir/00-basecall/processed", mode: 'copy', pattern: '*{fastq.gz,fastq,bam}'

    input:
    tuple val(samplename), val(barcode), path(bam_pass)
    
    output: 
    path('*{fastq.gz,fastq,bam}')
    
    script:
    """
    samtools cat ${bam_pass}/${barcode}/*.bam > ${samplename}.bam
    """
}

process READ_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/00-basecall", mode: 'copy', pattern: '*readstats.tsv'
    tag "${reads.simpleName}, ${reads.extension} file"

    input:
        path(reads)

    output:
        path("*readstats.tsv")

    script:
    
    """
    echo "file\treads\tbases\tn_bases\tmin_len\tmax_len\tn50\tGC_percent\tQ20_percent" > ${reads.simpleName}.readstats.tsv
    
    if [[ ${reads.extension} == bam ]]; then
        samtools fastq ${reads} | faster2 -ts - >> ${reads.simpleName}.readstats.tsv
    else 
    faster2 -ts ${reads} >> ${reads.simpleName}.readstats.tsv
    fi
    """
}
