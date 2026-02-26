def logColors() {
    return [
        reset: "\033[0m",
        green: "\033[0;32m",
        blue: "\033[0;34m",
        yellow: "\033[0;33m",
        cyan: "\033[0;36m"
    ]
}

def showHelp() {
    def c = logColors()

    log.info """
=============================================
${c.yellow}NXF-ALIGNMENT${c.reset}
${c.yellow}basecal, align, and analyze ONT data${c.reset}
=============================================

${c.yellow}Required/important options:${c.reset}
    ${c.green}--ref <path>${c.reset}           Reference FASTA (required unless --basecall or --report is used)

${c.yellow}Input options:${c.reset}
    ${c.green}--pod5 <dir>${c.reset}           Directory with POD5 files (use when basecalling)
    ${c.green}--samplename <str>${c.reset}     Sample name to use (if not provided, sample name is obtained from the pod5 file)
    ${c.green}--reads <file|dir>${c.reset}     BAM/FASTQ file or directory of reads (skips basecalling)
    ${c.green}--asfile <file>${c.reset}        Adaptive sampling CSV (filters reads to basecall)
    ${c.green}--kit${c.reset}                  Use for barcoded run - barcoding kit name (--samplesheet required)
    ${c.green}--samplesheet${c.reset}          Use for barcoded run - CSV or XLSX with columns: sample,barcode (--kit required)

${c.yellow}Output & config:${c.reset}
    ${c.green}-profile${c.reset}               Nextflow profile (standard, test, dev, singularity)
    ${c.green}--outdir${c.reset}               Output directory name (default: results)
    ${c.green}--basecall${c.reset}             Run the pipeline up to basecalling only
    ${c.green}--report${c.reset}               Run the pipeline up to reporting only (skips alignment and variants)
    ${c.green}--cpus${c.reset}                 Number of CPUs to use (default: 64)

${c.yellow}Processing options:${c.reset}
    ${c.green}--model${c.reset}                Dorado basecalling model (default: fast). For modifications use for example 'hac,5mCG_5hmCG'
    ${c.green}--herro${c.reset}                Enable herro correction (default: false). The corrected reads will be in 00-basecall, but will NOT be used in alignment.
    ${c.green}--bed${c.reset}                  BED file with regions (auto-generated from reference if omitted)
    ${c.green}--snp${c.reset}                  Enable SNP/small INDEL variant calling
    ${c.green}--snp_caller${c.reset}           SNP variant caller to use, only when --snp is specified (default: clair3, options: clair3, deepvariant)
    ${c.green}--deepvariant_model${c.reset}    DeepVariant model to use, only when --snp and --snp_caller deepvariant is specified (default: ONT_R104)
    ${c.green}--clair3_platform${c.reset}      Clair3 platform to use, only when --snp and --snp_caller clair3 is specified (default: ONT)
    ${c.green}--clair3_model${c.reset}         Clair3 model to use, only when --snp and --snp_caller clair3 is specified (default: r1041_e82_400bps_hac_v500)
    ${c.green}--sv${c.reset}                   Enable structural variant calling with sniffles2
    ${c.green}--phase${c.reset}                Enable SNP phasing with Whatshap (use only with --snp, only diploid cases supported)
    ${c.green}--annotate${c.reset}             Enable SNP variant annotation with snpEff (use only with --snp)
    ${c.green}--anno_db${c.reset}              snpEff database to use, only when --annotate is specified (default: hg38)
    ${c.green}--anno_filterQ${c.reset}         Filter out variants with quality lower than this before annotation (default: 20)
    ${c.green}--mods${c.reset}                 Perform base modifications analysis using modkit (requires --ref, read data must have mods)
    ${c.green}--mods_filter${c.reset}          Minimum coverage for base modifications calls (default: 5)
    """.stripIndent()
}