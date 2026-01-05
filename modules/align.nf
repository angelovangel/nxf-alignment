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

process BEDTOOLS_COMPLEMENT {
    container 'docker.io/biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    input:
        path(bed) 
        path(genome)

    output:
        path "*complement.bed"

    script:
    """
    sort -k1,1 -k2,2n ${bed} > ${bed.simpleName}.sorted.bed
    bedtools complement -i ${bed.simpleName}.sorted.bed -g ${genome} > ${bed.simpleName}.complement.bed
    """
}

process SAMTOOLS_BEDCOV {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${bam.simpleName}"

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*.bedcov.tsv'

    input:
        tuple path(bam), path(bai), path(bed), path(bedcomplement)

    output:
        path "*bedcov.tsv", emit: ch_bedcov
        path "*bedcov.compl.tsv", emit: ch_bedcov_complement
        path "*flagstat.json", emit: ch_flagstat

    script:
    """
    samtools bedcov ${bed} ${bam} > ${bam.simpleName}.bedcov.tsv
    samtools bedcov ${bedcomplement} ${bam} > ${bam.simpleName}.bedcov.compl.tsv
    samtools flagstat -O json ${bam} > ${bam.simpleName}.flagstat.json
    """
}

process REF_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${ref.simpleName}"

    input:
        path ref

    output:
        path "ref_stats.csv", emit: ch_ref_stats
        path "ref.genome", emit: ch_genome

    script:
    """
    samtools faidx ${ref}
    awk 'BEGIN {sum=0; count=0; print "contigs,bases"} {sum+=\$2; count++} END {print count "," sum}' ${ref}.fai > ref_stats.csv
    sort -k1,1 -k2,2n ${ref}.fai > ref.genome
    #cp ${ref}.fai ref.genome
    """
}