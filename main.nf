include {DORADO_BASECALL; DORADO_BASECALL_BARCODING;DORADO_CORRECT} from './modules/basecall.nf'
include {DORADO_ALIGN; MAKE_BEDFILE; BEDTOOLS_COV; BEDTOOLS_COMPLEMENT; SAMTOOLS_BEDCOV; DEEPTOOLS_BIGWIG; REF_STATS} from './modules/align.nf'
include {VCF_CLAIR3; VCF_DEEPVARIANT; VCF_STATS as VCF_STATS_SNP; VCF_STATS as VCF_STATS_SV; VCF_SNIFFLES2; VCF_PHASE; VCF_ANNOTATE; VCF_ANNOTATE_REPORT} from './modules/variants.nf'
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
    Profile             : ${workflow.profile}
    Container Engine    : ${workflow.containerEngine ?: 'local'}
    -------------------------------
    """.stripIndent()
    params.each { name, value -> log.info "${name.padRight(20)}: ${value}"}
log.info "=============================="

def showHelp() {
        log.info """
=============================================
NXF-ALIGNMENT
basecal, align, and analyze ONT data
=============================================

Required/important options:
    --ref <path>           Reference FASTA (required unless --basecall or --report is used)

Input options:
    --pod5 <dir>           Directory with POD5 files (use when basecalling)
    --reads <file|dir>     BAM/FASTQ file or directory of reads (skips basecalling)
    --asfile <file>        Adaptive sampling CSV (filters reads to basecall)

Processing options:
    --model                Dorado basecalling model (default: fast). For modifications use for example 'hac,5mCG_5hmCG'
    --herro                Enable herro correction (default: false). The corrected reads will be in 00-basecall, but will NOT be used in alignment.
    --kit                  Barcoding kit name (required with --samplesheet)
    --samplesheet          CSV with columns: sample,barcode (required with --kit)
    --bed                  BED file with regions (auto-generated from reference if omitted)
    --snp                  Enable SNP/small INDEL variant calling
    --snp_caller           SNP variant caller to use, only when --snp is specified (default: clair3, options: clair3, deepvariant)
    --deepvariant_model    DeepVariant model to use, only when --snp and --snp_caller deepvariant is specified (default: ONT_R104)
    --clair3_platform      Clair3 platform to use, only when --snp and --snp_caller clair3 is specified (default: ONT)
    --clair3_model         Clair3 model to use, only when --snp and --snp_caller clair3 is specified (default: r1041_e82_400bps_hac_v500)
    --sv                   Enable structural variant calling with sniffles2
    --phase                Enable SNP phasing with Whatshap (use only with --snp, only diploid cases supported)
    --annotate             Enable SNP variant annotation with snpEff (use only with --snp)
    --anno_db              snpEff database to use, only when --annotate is specified (default: hg38)
    --anno_filterQ         Filter out variants with quality lower than this before annotation (default: 20)

Output & config:
    --outdir               Output directory name (default: results)
    -profile               Nextflow profile (standard, test, dev, singularity)
    --basecall             Run the pipeline up to basecalling only
    --report               Run the pipeline up to reporting only (skips alignment and variants)

""".stripIndent()
}

// create empty placeholder files if not exist
/*
["runinfo", "refstats", "hist", "bedcov", "bedcov_compl", "flagstat", "variants", "sv_variants", "phase_stats"].each { name ->
    def f = file("${workflow.workDir}/empty_${name}" + (name.contains("info") || name.contains("stats") ? ".csv" : ""))
    if (!f.exists()) f.text = ""
    // assign to variables for easy reference
    if (name == "runinfo") empty_runinfo = f
    if (name == "refstats") empty_refstats = f
}
*/

// Additional specific assignments if needed
def empty_runinfo = file("${workflow.workDir}/empty_runinfo.csv")
def empty_refstats = file("${workflow.workDir}/empty_refstats.csv")
def empty_hist = file("${workflow.workDir}/empty_hist")
def empty_bedcov = file("${workflow.workDir}/empty_bedcov")
def empty_bedcov_compl = file("${workflow.workDir}/empty_bedcov_compl")
def empty_flagstat = file("${workflow.workDir}/empty_flagstat")
def empty_variants = file("${workflow.workDir}/empty_variants")
def empty_sv_variants = file("${workflow.workDir}/empty_sv_variants")
def empty_phase_stats = file("${workflow.workDir}/empty_phase_stats")

// Workflow properties - create CSV content as a string
def as_status = params.asfile ? "Yes" : "No"
def git_commit = "Unknown"
try {
    def git_proc = "git -C ${workflow.projectDir} rev-parse HEAD".execute()
    git_proc.waitFor()
    if (git_proc.exitValue() == 0) {
        git_commit = git_proc.text.trim().take(7)
    }
} catch (Exception e) {
    log.warn "Failed to get git commit: ${e.message}"
}
def workflow_properties = """\
CommandLine,User Name,Revision ID,Session ID,Run Name,Adaptive Sampling,Git Commit
"${workflow.commandLine}","${workflow.userName}","${workflow.scriptId.take(10)}","${workflow.sessionId.toString()}","${workflow.runName}","${as_status}","${git_commit}"
""".stripIndent()

// Create channel with CSV file
def ch_wf_properties = Channel.of(workflow_properties)
    .collectFile(name: 'wf_properties.csv', newLine: false)

// if no asfile, use dummy placeholder to still do dorado basecalling without as filtering
def ch_asfile = params.asfile ? Channel.fromPath(params.asfile, checkIfExists: true) : Channel.fromPath('EMPTY', type: 'file')

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

    if (params.kit) {
        DORADO_BASECALL_BARCODING(ch_asfile, ch_pod5)  
        
        ch_samplesheet
        .splitCsv(header:true)
        .filter{ it -> it.barcode =~ /^barcode*/ }
        .map { row -> tuple( row.sample, row.barcode ) }
        .combine( DORADO_BASECALL_BARCODING.out.ch_bam_pass )
        | MERGE_READS 

        if (params.herro) {
            DORADO_CORRECT(MERGE_READS.out)
        }
    } else {
        DORADO_BASECALL(ch_asfile, ch_pod5)
        // as an extra, do herro correction, the corrected reads are not used downstream
        if (params.herro) {
            DORADO_CORRECT(DORADO_BASECALL.out)
        }
    }
    
    emit: 
    ch_bc = params.kit ? MERGE_READS.out : DORADO_BASECALL.out
}

workflow {
    if (!params.basecall && !params.report && !params.ref) {
        error "Reference FASTA file (--ref) is required for full pipeline analysis (alignment and variant calling). Use --basecall or --report to run subsets of the pipeline without a reference."
    }

    // If 'reads' parameter is provided create a channel from that path.
    // also possible to pass a folder with reads, every read is one sample
    if (params.reads) {        
        if ( file(params.reads).isDirectory() ) {
            pattern = "*.{bam,fasta,fastq,fastq.gz,fq,fq.gz}"
            ch_reads = Channel.fromPath(params.reads + "/" + pattern, type: 'file', checkIfExists: true)           
        } else {
            ch_reads = Channel.fromPath(params.reads, checkIfExists: true)        
        }
    } else {
        // Otherwise, source the channel from the 'basecall' (or basecall + merge_reads) workflow's output.
        ch_reads = basecall().ch_bc
    }

    if (params.basecall) {
        return
    }

    if (params.ref) {
        ch_ref = Channel.fromPath(params.ref, checkIfExists: true)
        REF_STATS(ch_ref)
        ch_ref_stats = REF_STATS.out.ch_ref_stats
        ch_genome = REF_STATS.out.ch_genome
    } else {
        ch_ref = Channel.empty()
        ch_ref_stats = Channel.fromPath(empty_refstats)
        ch_genome = Channel.empty()
    }

    RUN_INFO( ch_reads.filter{ it.name.endsWith('.bam') }.first() )
    READ_STATS(ch_reads)

    if (params.report) {
        REPORT(
            RUN_INFO.out.ifEmpty(empty_runinfo),
            ch_wf_properties,
            READ_STATS.out.collect(),
            ch_ref_stats,
            //
            Channel.fromPath(empty_hist),
            Channel.fromPath(empty_bedcov),
            Channel.fromPath(empty_bedcov_compl),
            Channel.fromPath(empty_flagstat),
            Channel.fromPath(empty_variants),
            Channel.fromPath(empty_sv_variants),
            Channel.fromPath(empty_phase_stats),
            ch_asfile
        )
        return
    }

    // if no bedfile provided, just use the ref to generate one with the fasta entries
    if ( !params.bed ) {
        // generating bedfile from reference
        MAKE_BEDFILE(ch_ref)
        ch_bedfile = MAKE_BEDFILE.out
    } else {
        ch_bedfile = Channel.fromPath(params.bed, checkIfExists: true)
    }

    ch_ref \
    .combine( ch_reads ) \
    //.view()
    | DORADO_ALIGN \
    | combine( ch_bedfile ) \
    | BEDTOOLS_COV \
    
    DORADO_ALIGN.out
    .combine( ch_bedfile )
    .combine( BEDTOOLS_COMPLEMENT(ch_bedfile, ch_genome) )
    | SAMTOOLS_BEDCOV

    DORADO_ALIGN.out
    | DEEPTOOLS_BIGWIG

    ch_vc_input = DORADO_ALIGN.out
    .combine( ch_ref )
    .combine( ch_bedfile )

    // Variant Calling Logic
    if (params.sv) {
        ch_vc_input | VCF_SNIFFLES2
        VCF_STATS_SV(VCF_SNIFFLES2.out, 'sv')
    }
    // annotate and phase are within snp
    ch_vcf = Channel.empty()
    if (params.snp) {
        
        if (params.snp_caller == 'deepvariant') {
             ch_vc_input | VCF_DEEPVARIANT
             ch_vcf = VCF_DEEPVARIANT.out
        } else {
             ch_vc_input | VCF_CLAIR3
             ch_vcf = VCF_CLAIR3.out
        }
        
        VCF_STATS_SNP(ch_vcf, 'snp')

        if (params.annotate) {
            ch_vcf | VCF_ANNOTATE
            VCF_ANNOTATE.out.ch_vcfann_stats.collect() | VCF_ANNOTATE_REPORT
        }

        if (params.phase) {
            ch_bam_phase = DORADO_ALIGN.out.map{ bam, bai -> 
                def sample = bam.simpleName.replace('.align', '')
                tuple(sample, bam, bai)
            }

            ch_vcf_phase = ch_vcf.map{ vcf, tbi -> 
                def sample = vcf.simpleName.replace('.align.snp', '').replace('.snp', '')
                tuple(sample, vcf, tbi)
            }
            VCF_PHASE(
                ch_bam_phase.join(ch_vcf_phase),
                ch_ref.first(),
                ch_genome.first()
            )
            ch_phase_stats = VCF_PHASE.out.ch_vcfphase_stats.collect()
        } else {
            ch_phase_stats = Channel.fromPath(empty_phase_stats)
        }
    }
    
    REPORT(
        RUN_INFO.out.ifEmpty(empty_runinfo),
        ch_wf_properties,
        READ_STATS.out.collect(),
        ch_ref_stats,
        BEDTOOLS_COV.out.ch_hist.collect(),
        SAMTOOLS_BEDCOV.out.ch_bedcov.collect(),
        SAMTOOLS_BEDCOV.out.ch_bedcov_complement.collect(),
        SAMTOOLS_BEDCOV.out.ch_flagstat.collect(),
        params.snp ? VCF_STATS_SNP.out.collect() : Channel.fromPath(empty_variants),
        params.sv ? VCF_STATS_SV.out.collect() : Channel.fromPath(empty_sv_variants),
        params.phase ? ch_phase_stats : Channel.fromPath(empty_phase_stats),
        ch_asfile
    )  
}
