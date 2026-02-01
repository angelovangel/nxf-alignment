#!/usr/bin/env nextflow


process DORADO_BASECALL {
    // pod5 view is available in container
    container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/00-basecall", mode: 'copy'
    tag "${params.model}, ${params.asfile ? 'with' : 'no'} adaptive sampling"
    input:
        path decisionfile
        path pod5

    output:
        path "*.bam"

    script:
    def readids = params.asfile ? "--read-ids accepted_reads.txt" : ""
    def samplename = params.samplename ? "${params.samplename}" : ""
    """
    if [ "$decisionfile" != "EMPTY" ]; then
        echo "asfile provided, proceeding to basecall filtered reads."
        awk -F',' '\$2 == "sequence"' '$decisionfile' | cut -f1 -d, > accepted_reads.txt
        nreads_accept=\$(wc -l < accepted_reads.txt)
        nreads_total=\$(wc -l < '$decisionfile')
        echo -e "found \$nreads_accept out of \$nreads_total reads to basecall " 
    fi

    dorado basecaller $readids ${params.model} ${pod5} > reads.bam

    # check if reads.bam has reads and exit if no
    nreads=\$(wc -l < reads.bam)

    if [ "\$nreads" -eq 0 ]; then
        echo "No reads found in reads.bam, exiting." >&2
        exit 1
    fi
    
    samplename=\$(pod5 view --include "sample_id" ${pod5} | head -2 | tail -1)
    clean_name=\$(echo "\$samplename" | LC_ALL=C tr -dc '[:graph:]')

    if [ "${params.samplename}" != null ]; then
        clean_name=${params.samplename}
    fi
    echo "samplename: \$samplename"
    echo "clean name: \$clean_name"
    mv reads.bam \$clean_name.bam

    """
}

process DORADO_BASECALL_BARCODING {

    container 'docker.io/nanoporetech/dorado:latest'

    //publishDir "${params.outdir}/00-basecall", mode: 'copy'
    tag "${params.model}, ${params.asfile ? 'with' : 'no'} adaptive sampling"

    input:
        path decisionfile
        path pod5

    output:
        path "bam_pass", type: 'dir', emit: ch_bam_pass

    script:
    def readids = params.asfile ? "--read-ids accepted_reads.txt" : ""
    """
    if [ "$decisionfile" != "EMPTY" ]; then
        echo "asfile provided, proceeding to basecall filtered reads."
        awk -F',' '\$2 == "sequence"' '$decisionfile' | cut -f1 -d, > accepted_reads.txt
        nreads_accept=\$(wc -l < accepted_reads.txt)
        nreads_total=\$(wc -l < '$decisionfile')
        echo -e "Found \$nreads_accept out of \$nreads_total reads to basecall " 
    fi

    dorado basecaller $readids --kit-name ${params.kit} -o basecall-${params.model} ${params.model} ${pod5}
    # the folder with barcodes is basecall-sup/folder1/folder2/folder3/bam_pass
    [ -d "basecall-${params.model}" ] || { echo "Basecalling output folder empty!" >&2; exit 1; }
    ln -s basecall-${params.model}/*/*/*/bam_pass bam_pass
    """
}

process DORADO_CORRECT {
    container 'docker.io/nanoporetech/dorado:latest'
    tag "${bam.simpleName}"

    publishDir "${params.outdir}/00-basecall", mode: 'copy'

    input:
        path bam

    output:
        path "*.corr.fasta", emit: ch_corr_fasta
        path "*.log"

    script:
    """
    samtools fastq ${bam} > reads.fastq
    dorado correct reads.fastq > ${bam.baseName}.corr.fasta 2> ${bam.baseName}.dorado-correct.log
    """
}
