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
    printf "file\\treads\\tbases\\tn_bases\\tmin_len\\tmax_len\\tn50\\tGC_percent\\tQ20_percent\\tmods\\n" > ${reads.simpleName}.readstats.tsv
    
    mods="-"
    if [[ ${reads.extension} == "bam" ]]; then
        # Detect modifications
        bamtags=\$(samtools view ${reads} | head -n 1000 | grep -o "MM:Z:[^[:space:]]*" | cut -d',' -f1 | cut -d':' -f3 | sort -u)
        mods=\$(echo \$bamtags | parse_mods.py)

        samtools fastq ${reads} | faster2 -ts - | tr -d '\\n' >> ${reads.simpleName}.readstats.tsv
        echo -e "\\t\$mods" >> ${reads.simpleName}.readstats.tsv
    else 
        faster2 -ts ${reads} | tr -d '\\n' >> ${reads.simpleName}.readstats.tsv
        echo -e "\\t-" >> ${reads.simpleName}.readstats.tsv
    fi
    """
}
