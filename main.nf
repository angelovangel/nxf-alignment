include {DORADO_BASECALL; DORADO_BASECALL_BARCODING} from './modules/basecall.nf'
include {DORADO_ALIGN; MAKE_BEDFILE; BEDTOOLS_COV} from './modules/align.nf'
include {MERGE_READS; READ_STATS} from './modules/reads.nf'
include {RUN_INFO} from './modules/runinfo.nf'


if (params.help) {
        showHelp()
        exit 0
}

// 
def showHelp() {
        log.info """
=============================================
NXF-ALIGNMENT
basecal, align, and analyze ONT data
=============================================

Required/important options:
    --ref <path>     Reference FASTA (required unless -entry basecall is used)

Input options:
    --pod5 <dir>           Directory with POD5 files (use when basecalling)
    --reads <file|dir>     BAM/FASTQ file or directory of reads (skips basecalling)
    --asfile <file>        Adaptive sampling CSV (filters reads to basecall)

Processing options:
    --model <name>         Dorado basecalling model (default: fast)
    --kit <name>           Barcoding kit name (required with --samplesheet)
    --samplesheet <file>   CSV with columns: sample,barcode (required with --kit)
    --bed <file>           BED file with regions (auto-generated from reference if omitted)

Output & config:
    --outdir <name>        Output directory name (default: results)
    -profile <name>        Nextflow profile (test)
    -entry <name>          Workflow entry point (basecall, use for basecalling only)

""".stripIndent()
}


if (params.kit && !params.samplesheet) {
    error "If --kit is specified, --samplesheet must also be provided."
}

if (!params.kit && params.samplesheet) {
    error "If --samplesheet is provided, --kit must also be specified."
}





process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*html'
    
    input:
       tuple path(hist), path(readstats), path(runinfo)
    output:
        path "*.html"

    script:
    """
    make-report.py --hist $hist --readstats $readstats --runinfo $runinfo -o nxf-alignment-report.html
    """
}


workflow basecall {
    ch_pod5 = Channel.fromPath(params.pod5, checkIfExists: true)
    ch_samplesheet = params.samplesheet ? Channel.fromPath(params.samplesheet, checkIfExists: true) : null

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
    ch_ref = Channel.fromPath(params.ref)

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
    if ( !params.bed ) {
        // generating bedfile from reference
        MAKE_BEDFILE(Channel.fromPath(params.ref, checkIfExists: true))
        ch_bedfile = MAKE_BEDFILE.out
    } else {
        ch_bedfile = Channel.fromPath(params.bed, checkIfExists: true)
    }

    ch_reads
    .first() \
    | RUN_INFO

    ch_reads
    | READ_STATS

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
    .combine( RUN_INFO.out )
    //.view()
    | REPORT
    
}
