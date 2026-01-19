#!/usr/bin/env nextflow

process VCF_CLAIR3 {

    //container 'docker.io/hkubal/clair3:latest'
    container 'docker.io/hkubal/clair3-gpu:latest'

    publishDir "${params.outdir}/03-variants", mode: 'copy'
    errorStrategy 'ignore'
    tag "${bam.simpleName}, ${task.cpus} cpus"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.snp.vcf.gz"), path("${bam.simpleName}.snp.vcf.gz.tbi")

    script:
    def model = "${params.clair3_model}"
    def platform = "${params.clair3_platform}"
    """
    samtools faidx $ref

    /opt/bin/run_clair3.sh \
    --bam_fn=$bam \
    --ref_fn=$ref \
    --bed_fn=$bedfile \
    --platform=$platform \
    --model_path="/opt/models/$model" \
    --threads=${task.cpus} \
    --output="clair3_output"

    mv clair3_output/merge_output.vcf.gz ${bam.simpleName}.snp.vcf.gz
    mv clair3_output/merge_output.vcf.gz.tbi ${bam.simpleName}.snp.vcf.gz.tbi
    """ 
}

process VCF_DEEPVARIANT {
    //container 'docker.io/google/deepvariant:1.10.0-beta-gpu'
    publishDir "${params.outdir}/03-variants", mode: 'copy'
    tag "${bam.simpleName} ${task.cpus} cpus"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.snp.vcf.gz"), path("${bam.simpleName}.snp.vcf.gz.tbi")

    script:
    """
    samtools faidx $ref

    /opt/deepvariant/bin/run_deepvariant \
    --model_type=${params.deepvariant_model} \
    --ref=$ref \
    --reads=$bam \
    --regions=$bedfile \
    --output_vcf=${bam.simpleName}.snp.vcf.gz \
    --num_shards=${task.cpus}
    """ 
}

process VCF_STATS {
    
    container 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'
    errorStrategy 'ignore'
    tag "${vcf.simpleName}"

    input:
    tuple path(vcf), path(vcf_tbi)
    val vcftype

    output:
    path("${vcf.simpleName}.${vcftype}.query")

    script:
    def snp_header = "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%FILTER\t%QUAL\n"
    def sv_header = "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%FILTER\t%QUAL\t%INFO/SVTYPE\n"
    def format = vcftype == 'sv' ? sv_header : snp_header
    """
    bcftools query -f '$format' $vcf > ${vcf.simpleName}.${vcftype}.query
    """
}

process VCF_SNIFFLES2 {
    container 'docker.io/hydragenetics/sniffles2:2.6.3'
    publishDir "${params.outdir}/03-variants", mode: 'copy'
    tag "${bam.simpleName}"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.sv.vcf.gz"), path("${bam.simpleName}.sv.vcf.gz.tbi")

    script:
    """
    sniffles \
    --input $bam \
    --regions $bedfile \
    --vcf ${bam.simpleName}.sv.vcf.gz \
    --threads 4
    """ 
}

process VCF_ANNOTATE {
    container 'docker.io/dceoy/snpeff:latest'
    containerOptions '--entrypoint ""' 
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/04-annotations", mode: 'copy'
    tag "${vcf.simpleName} (filterQ >= ${params.anno_filterQ})"

    input:
    tuple path(vcf), path(vcf_tbi)

    output:
    path("${vcf.simpleName}.ann.vcf")
    path("${vcf.simpleName}.stats.csv"), emit: ch_vcfann_stats

    script:
    """
    # download snpEff database here
    mkdir -p ./snpeff_data

    SnpSift filter "(QUAL>=${params.anno_filterQ})" $vcf > filtered.vcf
    snpEff ann -dataDir \$PWD/snpeff_data -csvStats ${vcf.simpleName}.stats.csv ${params.anno_db} filtered.vcf > ${vcf.simpleName}.ann.vcf
    
    """
}

process VCF_ANNOTATE_REPORT {
    publishDir "${params.outdir}/05-annotations", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path(ann_stats)

    output:
    path("variants_annotation_report.html")

    script:
    """
    make-variants-report.py $ann_stats -o variants_annotation_report.html --filterQ ${params.anno_filterQ}
    """
}

process VCF_PHASE {
    container 'docker.io/tianjie16/whatshap:2.8'
    publishDir "${params.outdir}/04-phasing", mode: 'copy', pattern: "*.{vcf.gz,gtf,ht.bam,ht.bam.bai}"
    tag "${vcf.simpleName}"

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), path(vcf_tbi)
    path ref
    path genome

    output:
    tuple val(sample), path("${sample}.phase.vcf.gz"), path("${sample}.phase.gtf")
    path("${sample}.phase.tsv"), emit: ch_vcfphase_stats
    tuple path("${sample}.ht.bam"), path("${sample}.ht.bam.bai")

    script:
    """
    mv $genome ${ref}.fai

    whatshap phase \
    --reference $ref \
    --ignore-read-groups \
    -o ${sample}.phase.vcf.gz \
    $vcf \
    $bam 

    whatshap stats \
    --gtf ${sample}.phase.gtf \
    --tsv ${sample}.phase.tsv \
    ${sample}.phase.vcf.gz

    tabix ${sample}.phase.vcf.gz

    whatshap haplotag \
    -o ${sample}.ht.bam \
    --reference $ref \
    --ignore-read-groups \
    --output-threads=4 \
    ${sample}.phase.vcf.gz $bam

    samtools index ${sample}.ht.bam
    """
}
    