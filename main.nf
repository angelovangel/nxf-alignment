/*
* PARAMS
*/

params.pod5 = "${projectDir}/data/pod5"
params.asfile = "${projectDir}/data/as_decisions.csv"
params.model = "fast"

params.reference = "${projectDir}/data/ref.fasta"
params.bedfile = "${projectDir}/data/regions.bed"
params.outdir = "results"

params.reads = null // Parameter for reads when -entry align is used

// BASECALL

process DORADO_BASECALL {

    //container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/00-basecall", mode: 'copy'

    input:
        path decisionfile
        path pod5

    output:
        path "reads.bam"

    script:
    """
    awk -F',' '\$2 == "sequence"' '$decisionfile' | cut -f1 -d, > accepted_reads.txt
    nreads_accept=\$(wc -l < accepted_reads.txt)
    nreads_total=\$(wc -l < '$decisionfile')
    echo -e "Found \$nreads_accept out of \$nreads_total reads to basecall " 

    dorado basecaller -l accepted_reads.txt ${params.model} ${pod5} > reads.bam

    # check if reads.bam has reads and exit if no
    nreads=\$(samtools view -c reads.bam)

    if [ "\$nreads" -eq 0 ]; then
        echo "No reads found in reads.bam, exiting." >&2
        exit 1
    fi

    """
}

// ALIGN
process DORADO_ALIGN {

    container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/01-align", mode: 'copy', pattern: '*{bam,bai}'

    input:
        path ref 
        path reads

    output:
        tuple path('align.bam'), path('align.bam.bai')

    script:
    """
    dorado aligner ${ref} ${reads} | samtools sort -o align.bam
    samtools index align.bam
    """
}

// STATS
process SAMTOOLS_COV {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}/03-coverage", mode: 'copy', pattern: '*tsv'

    input:
        path bed
        tuple path(bam), path(bai)

    output:
        path "coverage.tsv"

    script:
    """
    echo -e "chr\tstart\tend\tlabel\tbases\tcoverage" > coverage.tsv
    samtools bedcov ${bed} ${bam} | awk '{
        region_length = \$3 - \$2
        coverage = (\$5 / region_length)
        print \$0 "\t" coverage
    }' >> coverage.tsv
    """
}

ch_ref = Channel.fromPath(params.reference, checkIfExists: true)

workflow basecall {
    ch_pod5 = Channel.fromPath(params.pod5, checkIfExists: true)
    ch_decisionfile = Channel.fromPath(params.asfile, checkIfExists: true)
    
    main:
    DORADO_BASECALL(ch_decisionfile, ch_pod5)

    emit: 
    ch_bc = DORADO_BASECALL.out
}

workflow align {
    
    if (params.reads_bam) {
        // If 'reads_bam' parameter is provided (e.g., via -entry align ... --reads <path>),
        // create a channel from that path.
        ch_reads = Channel.fromPath(params.reads, checkIfExists: true)
    } else {
        // Otherwise, source the channel from the 'basecall' workflow's output.
        ch_reads = basecall().ch_bc
    }
    
    main:
    DORADO_ALIGN(ch_ref, ch_reads)

    emit:
    ch_align = DORADO_ALIGN.out
}

workflow {
    align_out = align()
    SAMTOOLS_COV(params.bedfile, align.out.ch_align)
}
