#!/usr/bin/env nextflow

process DORADO_ALIGN {
    container 'docker.io/nanoporetech/dorado:latest'

    publishDir "${params.outdir}/01-align", mode: 'copy', pattern: "*.align.bam*"
    tag "${sample}, ${task.cpus} cpus"

    input:
        tuple val(sample), path(ref), path(reads)

    output:
        tuple val(sample), path("${sample}.align.bam"), path("${sample}.align.bam.bai")
        path "versions.txt", emit: versions

    script:
    """
    # Calculate a balanced distribution of threads to avoid oversubscribing task.cpus
    FASTQ_CPUS=\$(( ${task.cpus} > 8 ? 4 : 1 ))
    SORT_CPUS=\$(( ${task.cpus} > 4 ? ${task.cpus} / 4 : 1 ))
    MAP_CPUS=\$(( ${task.cpus} - \$FASTQ_CPUS - \$SORT_CPUS ))
    [[ \$MAP_CPUS -lt 1 ]] && MAP_CPUS=1

    samtools fastq -@ \$FASTQ_CPUS -T '*' ${reads} | \
    minimap2 ${params.mmap2} -t \$MAP_CPUS -y ${ref} - | \
    samtools sort -@ \$SORT_CPUS -m 2G -T "sort_tmp" -o ${sample}.align.bam -

    # exit here if no reads map to ref
    nreads=\$(samtools view -@ ${task.cpus} -c ${sample}.align.bam)
    if [ "\$nreads" -eq 0 ]; then
        echo "No reads mapped to reference. Exiting."
        exit 1
    fi

    samtools index -@ ${task.cpus} ${sample}.align.bam

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: minimap2 v\$(minimap2 --version)
    ${task.process}: samtools v\$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
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

process BEDTOOLS_COV {
    container 'docker.io/biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*hist.tsv'
    publishDir "${params.outdir}/logs", mode: 'copy', pattern: '*.bed'
    tag "${sample}"

    input:
        tuple val(sample), path(bam), path(bai), path(bed)

    output:
        path "*hist.tsv", emit: ch_hist
        path "nxf-alignment.bed"
        path "versions.txt", emit: versions

    script:
    """
    echo -e "chr\tstart\tend\tlabel\tdepth\tbases_at_depth\tsize\tpercent_at_depth" > ${sample}.hist.tsv
    bedtools coverage -a ${bed} -b ${bam} -hist >> ${sample}.hist.tsv
    cp $bed nxf-alignment.bed

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: bedtools v\$(bedtools --version | sed 's/^bedtools v//')
    END_VERSIONS
    """
}

process BEDTOOLS_COMPLEMENT {
    container 'docker.io/biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

    input:
        path(bed) 
        path(genome)

    output:
        path "*complement.bed"

    script:
    """
    sort -k1,1 -k2,2n ${bed} > ${bed.simpleName}.sorted.bed
    bedtools complement -i ${bed.simpleName}.sorted.bed -g ${genome} > ${bed.simpleName}.complement.bed
    """
}

process SAMTOOLS_BEDCOV {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${sample}"

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*.bedcov.tsv'

    input:
        tuple val(sample), path(bam), path(bai), path(bed), path(bedcomplement)

    output:
        path "*bedcov.tsv", emit: ch_bedcov
        path "*bedcov.compl.tsv", emit: ch_bedcov_complement
        path "*flagstat.json", emit: ch_flagstat

    script:
    """
    samtools bedcov -c ${bed} ${bam} > ${sample}.bedcov.tsv
    samtools bedcov -c ${bedcomplement} ${bam} > ${sample}.bedcov.compl.tsv
    samtools flagstat -O json ${bam} > ${sample}.flagstat.json
    """
}

process DEEPTOOLS_BIGWIG {
    container 'mgibio/deeptools:3.5.3'
    tag "${sample}"

    publishDir "${params.outdir}/02-coverage", mode: 'copy', pattern: '*.bigwig'

    input:
        tuple val(sample), path(bam), path(bai)

    output:
        path "*.bigwig"
        path "versions.txt", emit: versions

    script:
    """
    bamCoverage -b ${bam} -o ${sample}.bigwig -p ${task.cpus}

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: deeptools v\$(bamCoverage --version | sed 's/^bamCoverage //')
    END_VERSIONS
    """
}

process REF_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${ref.simpleName}"

    input:
        path ref

    output:
        path "ref_stats.csv", emit: ch_ref_stats
        path "ref.genome", emit: ch_genome

    script:
    """
    samtools faidx ${ref}
    awk 'BEGIN {sum=0; count=0; print "contigs,bases"} {sum+=\$2; count++} END {print count "," sum}' ${ref}.fai > ref_stats.csv
    sort -k1,1 -k2,2n ${ref}.fai > ref.genome
    #cp ${ref}.fai ref.genome
    """
}