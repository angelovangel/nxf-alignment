process PGX_PANNO {
    container 'docker.io/aangeloo/panno:0.3.1'
    publishDir "${params.outdir}/04-pgx", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path(vcf), path(vcf_tbi)
    val population

    output:
    path("*.html")
    path "versions.txt", emit: versions

    script:
    """
    # PAnno requires decompressed VCF and crashes on missing genotypes (./.)
    if [[ "$vcf" == *.gz ]]; then
        zcat $vcf > unfiltered.vcf
    else
        cp $vcf unfiltered.vcf
    fi
    
    grep "^#" unfiltered.vcf > ${sample}.vcf
    grep -v "^#" unfiltered.vcf | awk '\$10 !~ /^[.]/ {print}' >> ${sample}.vcf

    panno -i ${sample}.vcf -p $population -s $sample -o .

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: panno \$(panno --version 2>&1 | head -n 1 | sed 's/^panno //')
    END_VERSIONS
    """
}

process PHARMCAT_PREPARE {
    container 'docker.io/pgkb/pharmcat:latest'
    
    output:
    path "pharmcat_positions.vcf.gz", emit: vcf
    path "pharmcat_positions.vcf.gz.tbi", emit: tbi

    script:
    """
    cp /pharmcat/pharmcat_positions.vcf.bgz pharmcat_positions.vcf.gz
    tabix pharmcat_positions.vcf.gz
    """
}

process PHARMCAT_GENOTYPE {
    container 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'
    tag "${sample}"

    input:
    tuple val(sample), path(bam), path(bai)
    path(ref)
    path(fai)
    path(positions_vcf)
    path(positions_tbi)

    output:
    tuple val(sample), path("${sample}.pharmcat_genotypes.vcf.gz"), path("${sample}.pharmcat_genotypes.vcf.gz.tbi"), emit: vcf
    path "versions.txt", emit: versions

    script:
    """
    # Create a BED file for available chromosomes from the reference index
    cut -f1 $fai > ref_chroms.txt
    awk -v OFS='\t' '{print \$1, 0, \$2}' $fai > ref_chroms.bed

    # Determine if we need to rename and/or filter the positions VCF
    if grep -q "^chr" ref_chroms.txt; then
        # Reference uses 'chr', so just filter positions by available chromosomes
        bcftools view -T ref_chroms.bed $positions_vcf -Oz -o local_positions.vcf.gz
    else
        # Reference does NOT use 'chr', strip 'chr' from positions and filter
        bcftools query -f '%CHROM\\n' $positions_vcf | sort -u | while read -r c; do
            echo "\$c \${c#chr}"
        done > chrom_map.txt
        bcftools annotate --rename-chrs chrom_map.txt $positions_vcf -Ou | \\
        bcftools view -T ref_chroms.bed -Oz -o local_positions.vcf.gz
    fi
    bcftools index -t local_positions.vcf.gz
    pos_vcf="local_positions.vcf.gz"

    set -o pipefail
    bcftools mpileup -f $ref -R \$pos_vcf $bam -a FORMAT/DP,FORMAT/AD | \\
    bcftools call -m -Ou | \\
    bcftools view -i 'FORMAT/DP>0' -Oz -o ${sample}.pharmcat_genotypes.vcf.gz
    bcftools index -t ${sample}.pharmcat_genotypes.vcf.gz

    echo "${task.process}: bcftools v\$(bcftools --version | head -n 1 | sed 's/^bcftools //')" > versions.txt
    """
}

process PHARMCAT_MERGE {
    container 'quay.io/biocontainers/bcftools:1.23--h3a4d415_0'
    tag "${sample}"

    input:
    tuple val(sample), path(variants_vcf), path(variants_tbi), path(genotypes_vcf), path(genotypes_tbi)
    path ref
    path fai

    output:
    tuple val(sample), path("${sample}.merged_pgx.vcf.gz"), path("${sample}.merged_pgx.vcf.gz.tbi"), emit: vcf
    path "versions.txt", emit: versions

    script:
    """
    # ensure sample names are correct
    echo "$sample" > sample.txt
    bcftools reheader -s sample.txt -o variants_rh.vcf.gz $variants_vcf
    bcftools reheader -s sample.txt -o genotypes_rh.vcf.gz $genotypes_vcf
    bcftools index -t variants_rh.vcf.gz
    bcftools index -t genotypes_rh.vcf.gz

    # concat prioritizing variants_rh (variants first)
    bcftools concat -a variants_rh.vcf.gz genotypes_rh.vcf.gz | \
    bcftools norm -c ws -f $ref -m -any -Ou | \
    bcftools sort -Ou | \
    bcftools norm -d all -Oz -o ${sample}.merged_pgx.vcf.gz
    bcftools index -t ${sample}.merged_pgx.vcf.gz

    echo "${task.process}: bcftools v\$(bcftools --version | head -n 1 | sed 's/^bcftools //')" > versions.txt
    """
}

process PGX_PHARMCAT {
    container 'docker.io/pgkb/pharmcat:latest'
    publishDir "${params.outdir}/04-pgx", mode: 'copy'
    tag "${sample}"

    input:
    tuple val(sample), path(vcf), path(vcf_tbi)
    path ref
    path fai

    output:
    path "${sample}_pharmcat/*"
    path "${sample}.pharmcat.html"
    path "${sample}.pgx-report.html"
    path "versions.txt", emit: versions

    script:
    """
    # Manually preprocess to ensure the VCF is sorted after normalization
    # PharmCAT's internal preprocessor sometimes fails to produce a sorted VCF after left-alignment
    bcftools view -s $sample -R /pharmcat/pharmcat_regions.bed $vcf | \\
    bcftools norm -f $ref -m -any -c ws -Ou | \\
    bcftools sort -Oz -o ${sample}.preprocessed.vcf.gz
    bcftools index -t ${sample}.preprocessed.vcf.gz

    # Run the core PharmCAT tool on the sorted, normalized VCF
    pharmcat \\
        -vcf ${sample}.preprocessed.vcf.gz \\
        -o ${sample}_pharmcat \\
        -bf ${sample} \\
        -matcher \\
        -phenotyper \\
        -reporter \\
        -reporterHtml \\
        -reporterJson
    
    if [[ -f "${sample}_pharmcat/${sample}.report.html" ]]; then
        mv ${sample}_pharmcat/${sample}.report.html ${sample}.pharmcat.html
        echo "PharmCAT report generated: ${sample}.pharmcat.html"
        # generate report
        vcf-report.py --vcf ${sample}.preprocessed.vcf.gz --report ${sample}_pharmcat/${sample}.report.json --out ${sample}.pgx-report.html
    else
        echo "PharmCAT report could not be generated for sample: ${sample}"
        exit 1
    fi

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: \$(pharmcat --version 2>&1 | head -n 1)
    END_VERSIONS
    """
}