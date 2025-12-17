include {DORADO_BASECALL; DORADO_BASECALL_BARCODING} from './modules/basecall.nf'
include {DORADO_ALIGN; MAKE_BEDFILE; BEDTOOLS_COV; BEDTOOLS_COMPLEMENT; SAMTOOLS_BEDCOV; REF_STATS} from './modules/align.nf'
include {MERGE_READS; READ_STATS} from './modules/reads.nf'
include {RUN_INFO} from './modules/runinfo.nf'
include {REPORT} from './modules/report.nf'

if (params.help) {
        showHelp()
        exit 0
}

// Log execution environment
log.info """
    NXF-ALIGNMENT Execution Summary
    ===============================
    Profile         : ${workflow.profile}
    Container Engine: ${workflow.containerEngine ?: 'local'}
    
    pod5            : ${params.pod5}
    reads           : ${params.reads}
    asfile          : ${params.asfile}
    model           : ${params.model}
    ref             : ${params.ref}
    bed             : ${params.bed}
    kit             : ${params.kit}
    samplesheet     : ${params.samplesheet}
    outdir          : ${params.outdir}
    ===============================
""".stripIndent()

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
    --model <name>         Dorado basecalling model (default: fast). For modifications use for example 'hac,5mCG_5hmCG'
    --kit <name>           Barcoding kit name (required with --samplesheet)
    --samplesheet <file>   CSV with columns: sample,barcode (required with --kit)
    --bed <file>           BED file with regions (auto-generated from reference if omitted)

Output & config:
    --outdir <name>        Output directory name (default: results)
    -profile <name>        Nextflow profile (test)
    -entry <name>          Workflow entry point (basecall - basecalling only, report - basecalling + report)

""".stripIndent()
}

// create empty runinfo file if not created
def empty_runinfo = file("${workflow.workDir}/empty_runinfo.csv")
empty_runinfo.text = ""

def empty_refstats = file("${workflow.workDir}/empty_refstats.csv")
empty_refstats.text = ""

// Workflow properties - create CSV content as a string
def workflow_properties = """\
CommandLine,UserName,RevisionID,SessionID,RunName
${workflow.commandLine},${workflow.userName},${workflow.scriptId.take(10)},${workflow.sessionId.toString()},${workflow.runName}
""".stripIndent()

// Create channel with CSV file
def ch_wf_properties = Channel.of(workflow_properties)
    .collectFile(name: 'wf_properties.csv', newLine: false)

if (params.kit && !params.samplesheet) {
    error "If --kit is specified, --samplesheet must also be provided."
}

if (!params.kit && params.samplesheet) {
    error "If --samplesheet is provided, --kit must also be specified."
}

// do basecall only
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
// do basecall + reporting
workflow report {
    // Calc ref stats if ref exists, else empty
    if (params.ref) {
        REF_STATS(Channel.fromPath(params.ref))
        ch_ref_stats = REF_STATS.out.ch_ref_stats
        
    } else {
        ch_ref_stats = Channel.fromPath(empty_refstats)
    }

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
    
    RUN_INFO( ch_reads.filter{ it.name.endsWith('.bam') }.first() )
    READ_STATS(ch_reads)
    
    REPORT(
        RUN_INFO.out.ifEmpty(empty_runinfo),
        ch_wf_properties,
        READ_STATS.out.collect(),
        ch_ref_stats,
        Channel.fromPath('./work/empty.hist.tsv', type: 'file'),
        Channel.fromPath('./work/empty.bedcov.tsv', type: 'file'),
        Channel.fromPath('./work/empty.bedcov_compl.tsv', type: 'file'),
        Channel.fromPath('./work/empty.flagstat.json', type: 'file'),
    )
}

workflow {
    ch_ref = Channel.fromPath(params.ref)
    REF_STATS(ch_ref)

    // If 'reads' parameter is provided create a channel from that path.
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

    RUN_INFO( ch_reads.filter{ it.name.endsWith('.bam') }.first() )
    READ_STATS(ch_reads)

    ch_ref \
    .combine( ch_reads ) \
    //.view()
    | DORADO_ALIGN \
    | combine( ch_bedfile ) \
    | BEDTOOLS_COV \
    
    DORADO_ALIGN.out
    .combine( ch_bedfile )
    .combine( BEDTOOLS_COMPLEMENT(ch_bedfile, REF_STATS.out.ch_genome) )
    | SAMTOOLS_BEDCOV

    REPORT(
        RUN_INFO.out.ifEmpty(empty_runinfo),
        ch_wf_properties,
        READ_STATS.out.collect(),
        REF_STATS.out.ch_ref_stats,
        BEDTOOLS_COV.out.collect(),
        SAMTOOLS_BEDCOV.out.ch_bedcov.collect(),
        SAMTOOLS_BEDCOV.out.ch_bedcov_complement.collect(),
        SAMTOOLS_BEDCOV.out.ch_flagstat.collect(),
    )  
}
