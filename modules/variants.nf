#!/usr/bin/env nextflow

process VCF_CLAIR3 {

    //container 'docker.io/hkubal/clair3:latest'
    container 'docker.io/hkubal/clair3-gpu:v1.2.0'

    publishDir "${params.outdir}/03-variants/snps", mode: 'copy', pattern: "*.{vcf.gz,vcf.gz.tbi}"
    errorStrategy 'ignore'
    tag "${bam.simpleName}, ${task.cpus} cpus"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.snp.vcf.gz"), path("${bam.simpleName}.snp.vcf.gz.tbi")
    path "versions.txt", emit: versions

    script:
    def model = "${params.clair3_model}"
    def platform = "${params.clair3_platform}"
    """
    samtools faidx $ref

    /opt/bin/run_clair3.sh \
    --bam_fn=$bam \
    --ref_fn=$ref \
    --bed_fn=$bedfile \
    --sample_name=${bam.simpleName} \
    --platform=$platform \
    --model_path="/opt/models/$model" \
    --threads=${task.cpus} \
    --output="clair3_output"

    mv clair3_output/merge_output.vcf.gz ${bam.simpleName}.snp.vcf.gz
    mv clair3_output/merge_output.vcf.gz.tbi ${bam.simpleName}.snp.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: clair3 v\$(/opt/bin/run_clair3.sh --version 2>&1 | sed 's/^Clair3 v//')
    END_VERSIONS
    """ 
}

process VCF_DEEPVARIANT {
    //container 'docker.io/google/deepvariant:1.10.0-beta-gpu'
    publishDir "${params.outdir}/03-variants/snps", mode: 'copy'
    tag "${bam.simpleName} ${task.cpus} cpus"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.snp.vcf.gz"), path("${bam.simpleName}.snp.vcf.gz.tbi")
    path "versions.txt", emit: versions

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

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: deepvariant v\$(/opt/deepvariant/bin/run_deepvariant --version 2>&1 | sed 's/^DeepVariant //')
    END_VERSIONS
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
    path "versions.txt", emit: versions

    script:
    def snp_header = "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%FILTER\t%QUAL\n"
    def sv_header = "%CHROM\t%POS\t%REF\t%ALT\t%TYPE\t%FILTER\t%QUAL\t%INFO/SVTYPE\n"
    def format = vcftype == 'sv' ? sv_header : snp_header
    """
    bcftools query -f '$format' $vcf > ${vcf.simpleName}.${vcftype}.query
    
    echo "${task.process}: bcftools v\$(bcftools --version | head -n 1 | sed 's/^bcftools //')" > versions.txt
    """
}

process VCF_SNIFFLES2 {
    container 'docker.io/hydragenetics/sniffles2:2.6.3'
    publishDir "${params.outdir}/03-variants/sv", mode: 'copy', pattern: "*.{vcf.gz,vcf.gz.tbi}"
    tag "${bam.simpleName}"

    input:
    tuple path(bam), path(bai), path(ref), path(bedfile)

    output:
    tuple path("${bam.simpleName}.sv.vcf.gz"), path("${bam.simpleName}.sv.vcf.gz.tbi")
    path "versions.txt", emit: versions

    script:
    """
    sniffles \
    --input $bam \
    --regions $bedfile \
    --vcf ${bam.simpleName}.sv.vcf.gz \
    --threads ${task.cpus}

    echo "${task.process}: sniffles2 v\$(sniffles --version | head -n 1 | sed 's/^Sniffles2, Version //')" > versions.txt
    """ 
}

process VCF_ANNOTATE {
    container 'docker.io/dceoy/snpeff:latest'
    containerOptions '--entrypoint ""' 
    errorStrategy 'ignore'
    
    publishDir "${params.outdir}/03-variants/annotations", mode: 'copy', pattern: "*.{ann.vcf,stats.csv}"
    tag "${vcf.simpleName} (filterQ >= ${params.anno_filterQ})"

    input:
    tuple path(vcf), path(vcf_tbi)

    output:
    path("${vcf.simpleName}.ann.vcf")
    path("${vcf.simpleName}.stats.csv"), emit: ch_vcfann_stats
    path "versions.txt", emit: versions

    script:
    """
    # download snpEff database here
    mkdir -p ./snpeff_data

    SnpSift filter "(QUAL>=${params.anno_filterQ})" $vcf > filtered.vcf
    snpEff ann -dataDir \$PWD/snpeff_data -csvStats ${vcf.simpleName}.stats.csv ${params.anno_db} filtered.vcf > ${vcf.simpleName}.ann.vcf

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: snpEff v\$(snpEff -version | head -n 1 | sed 's/^SnpEff //')
    END_VERSIONS
    """
}

process VCF_BGZIP {
    container 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'
    //publishDir "${params.outdir}/03-variants/annotations", mode: 'copy'
    tag "${vcf.simpleName}"

    input:
    path vcf

    output:
    tuple path("${vcf.simpleName}.vcf.gz"), path("${vcf.simpleName}.vcf.gz.tbi")

    script:
    """
    bgzip -c $vcf > ${vcf.simpleName}.vcf.gz
    tabix ${vcf.simpleName}.vcf.gz
    """
}

process MERGE_VARIANTS {
    container 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'
    publishDir "${params.outdir}/03-variants/merged", mode: 'copy'
    tag "${snp_vcf.simpleName}"

    input:
    tuple path(snp_vcf), path(snp_tbi), path(sv_vcf), path(sv_tbi)

    output:
    tuple path("*.merged.vcf.gz"), path("*.merged.vcf.gz.tbi")

    script:
    def sample = snp_vcf.simpleName.replace('.snp', '').replace('.align', '').replace('.ann', '').replace('.phase', '').replace('.vcf', '').replace('.gz', '')
    """
    echo "$sample" > sample_name.txt
    
    bcftools reheader -s sample_name.txt -o snp_reheaded.vcf.gz $snp_vcf
    bcftools reheader -s sample_name.txt -o sv_reheaded.vcf.gz $sv_vcf
    
    bcftools index snp_reheaded.vcf.gz
    bcftools index sv_reheaded.vcf.gz

    bcftools concat -a snp_reheaded.vcf.gz sv_reheaded.vcf.gz | bcftools sort -o ${sample}.merged.vcf.gz -O z
    tabix ${sample}.merged.vcf.gz
    """
}

process VCF_ANNOTATE_REPORT {
    publishDir "${params.outdir}", mode: 'copy'
    errorStrategy 'ignore'

    input:
    path(ann_stats)

    output:
    path("variants-annotation-report.html")

    script:
    """
    make-variants-report.py $ann_stats -o variants-annotation-report.html --filterQ ${params.anno_filterQ}
    """
}

process VCF_PHASE {
    container 'docker.io/tianjie16/whatshap:2.8'
    publishDir "${params.outdir}/03-variants/phasing", mode: 'copy', pattern: "*.{vcf.gz,vcf.gz.tbi,gtf}"
    publishDir "${params.outdir}/01-align", mode: 'copy', pattern: "*.{ht.bam,ht.bam.bai}"
    tag "${vcf.simpleName}"

    input:
    tuple val(sample), path(bam), path(bai), path(vcf), path(vcf_tbi)
    path ref
    path genome

    output:
    tuple val(sample), path("${sample}.phase.{vcf.gz,vcf.gz.tbi}"), path("${sample}.phase.gtf")
    path("${sample}.phase.tsv"), emit: ch_vcfphase_stats
    tuple path("${sample}.ht.bam"), path("${sample}.ht.bam.bai")
    path "versions.txt", emit: versions

    script:
    """
    mv $genome ${ref}.fai

    whatshap phase \
    --reference $ref \
    --ignore-read-groups \
    --tag HP \
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
    --output-threads ${task.cpus} \
    ${sample}.phase.vcf.gz $bam

    samtools index -@ ${task.cpus} ${sample}.ht.bam

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: whatshap v\$(whatshap --version)
    ${task.process}: samtools v\$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
    