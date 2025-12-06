


if (params.kit && !params.samplesheet) {
    error "If --kit is specified, --samplesheet must also be provided."
}

// BASECALL

process DORADO_BASECALL {

    //container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/00-basecall", mode: 'copy'
    tag "${params.asfile ? 'with' : 'no'} adaptive sampling"
    input:
        path decisionfile
        path pod5

    output:
        path "reads.bam"

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
    
    dorado basecaller $readids ${params.model} ${pod5} > reads.bam

    # check if reads.bam has reads and exit if no
    nreads=\$(wc -l < reads.bam)

    if [ "\$nreads" -eq 0 ]; then
        echo "No reads found in reads.bam, exiting." >&2
        exit 1
    fi

    """
}

process DORADO_BASECALL_BARCODING {

    //container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/00-basecall", mode: 'copy'
    tag "${params.asfile ? 'with' : 'no'} adaptive sampling"

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

process MERGE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    errorStrategy 'ignore' //because some barcodes defined in the samplesheet might be missing in the data
    
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
    //publishDir "${params.outdir}/00-basecall", mode: 'copy', pattern: '*readstats.tsv'
    tag "${reads.simpleName}"

    input:
        path(reads)

    output:
        path("*readstats.tsv")

    script:
    """
    echo "file\treads\tbases\tn_bases\tmin_len\tmax_len\tn50\tGC_percent\tQ20_percent" > ${reads.simpleName}.readstats.tsv
    samtools fastq ${reads} | faster2 -ts - >> ${reads.simpleName}.readstats.tsv
    """
}

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

/*
process SAMTOOLS_BEDCOV {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*cov.tsv'
    tag "${bam.simpleName}"

    input:
        tuple path(bam), path(bai), path(bed)

    output:
        path "*cov.tsv"

    script:
    """
    echo -e "chr\tstart\tend\tlabel\tbases\tregion_len\tcoverage" > ${bam.simpleName}.cov.tsv
    samtools bedcov ${bed} ${bam} | awk '{
        region_length = \$3 - \$2
        coverage = (\$5 / region_length)
        print \$0 "\t" region_length "\t" coverage
    }' >> ${bam.simpleName}.cov.tsv
    """
}
*/

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


process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*html'
    
    input:
       tuple path(hist), path(readstats)
    output:
        path "*.html"

    script:
    """
    
    make-report.py --hist $hist --readstats $readstats -o nxf-ontas-report.html
    """
}

ch_ref = Channel.fromPath(params.reference, checkIfExists: true)
ch_samplesheet = params.samplesheet ? Channel.fromPath(params.samplesheet, checkIfExists: true) : null

workflow basecall {
    ch_pod5 = Channel.fromPath(params.pod5, checkIfExists: true)
    // if no asfile, use dummy placeholder to still do dorado basecalling without as filtering
    ch_decisionfile = params.asfile ? Channel.fromPath(params.asfile, checkIfExists: true) : Channel.fromPath('EMPTY', type: 'file')
    
    if (params.kit) {
        DORADO_BASECALL_BARCODING(ch_decisionfile, ch_pod5)  
        
        ch_samplesheet
        .splitCsv(header:true)
        .filter{ it -> it.barcode =~ /^barcode*/ }
        .map { row -> tuple( row.sample, row.barcode ) }
        .combine( DORADO_BASECALL_BARCODING.out.ch_bam_pass )
        | MERGE_READS 
    } else {
        DORADO_BASECALL(ch_decisionfile, ch_pod5)
    }
    
    emit: 
    ch_bc = params.kit ? MERGE_READS.out : DORADO_BASECALL.out
}

workflow {
    
    // If 'reads' parameter is provided create a 
    // channel from that path.
    // also possible to pass a folder with reads, every read is one sample
    if (params.reads) {
        
        if ( file(params.reads).isDirectory() ) {
            pattern = "*.{bam,fastq,fastq.gz,fq,fq.gz}"
            ch_reads = Channel.fromPath(params.reads + "/" + pattern, type: 'file', checkIfExists: true)
            
        } else {
            ch_reads = Channel.fromPath(params.reads, checkIfExists: true)        
        }
    } else {
        // Otherwise, source the channel from the 'basecall' (or basecall + merge_reads) workflow's output.
        ch_reads = basecall().ch_bc
    }

    // if no bedfile provided, just use the ref to generate one with the fasta entries
    if ( !params.bedfile ) {
        // generating bedfile from reference
        MAKE_BEDFILE(Channel.fromPath(params.reference, checkIfExists: true))
        ch_bedfile = MAKE_BEDFILE.out
    } else {
        ch_bedfile = Channel.fromPath(params.bedfile, checkIfExists: true)
    }
    
    ch_reads
    | READ_STATS \
    
    ch_ref \
    .combine( ch_reads ) \
    //.view()
    | DORADO_ALIGN \
    | combine( ch_bedfile ) \
    | BEDTOOLS_COV \
    
    
    BEDTOOLS_COV.out
    .collect()
    .toList()
    .combine( READ_STATS.out.collect().toList() )
    //.view()
    | REPORT
    
}

