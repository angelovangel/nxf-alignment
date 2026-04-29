process VCF_PANNO {
    container 'docker.io/aangeloo/panno:0.3.1'
    publishDir "${params.outdir}/03-variants/pgx", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path(vcf), path(vcf_tbi)
    val population

    output:
    path("*.html")
    path "versions.txt", emit: versions

    script:
    """
    # PAnno requires decompressed VCF
    if [[ "$vcf" == *.gz ]]; then
        zcat $vcf > ${sample}.vcf
    else
        cp $vcf ${sample}.vcf
    fi

    panno -i ${sample}.vcf -p $population -s $sample -o .

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: panno \$(panno --version 2>&1 | head -n 1 | sed 's/^panno //')
    END_VERSIONS
    """
}

process VCF_PHARMCAT {
    container 'docker.io/pgkb/pharmcat:latest'
    publishDir "${params.outdir}/03-variants/pgx", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path(vcf), path(vcf_tbi)

    output:
    path "${sample}_pharmcat/*"
    path "${sample}.pharmcat.html"
    path "versions.txt", emit: versions

    script:
    """
    pharmcat_pipeline $vcf \
        -o ${sample}_pharmcat \
        -bf ${sample} \
        -reporterHtml \
        -reporterJson
    mv ${sample}_pharmcat/${sample}.report.html ${sample}.pharmcat.html

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: \$(pharmcat --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}