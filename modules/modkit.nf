
process MODKIT {
    container 'docker.io/ontresearch/modkit:latest'
    publishDir "${params.outdir}/04-modifications", mode: 'copy'
    errorStrategy { 
        if (task.exitStatus == 42) {
            println "MODKIT: [WARNING] No modified reads detected for ${task.tag}. Skipping."
            return 'ignore'
        }
        return 'ignore'
    }
    tag "${bam.simpleName}"

    input:
        tuple path(bam), path(bai)

    output:
    tuple path("*.bedmethyl.gz"), path("*.bedmethyl.gz.tbi")
    path "*.summary.tsv"
    path "modprobs"
    path "versions.txt", emit: versions

    script:
    def filter = params.mods_filter
    """
    cat <<-END_VERSIONS > versions.txt
    ${task.process}: modkit v\$(modkit --version 2>&1 | head -n 1 | sed 's/^modkit //')
    ${task.process}: tabix v\$(tabix --version 2>&1 | head -n 1 | sed 's/^tabix (htslib) //')
    END_VERSIONS

    modkit modbam summary ${bam} > ${bam.simpleName}.summary.tsv
    nreads=\$(grep "# total_reads_used" ${bam.simpleName}.summary.tsv | awk '{print \$NF}')

    if [ "\$nreads" -gt 0 ]; then
        modkit sample-probs -o modprobs --hist --prefix ${bam.simpleName} ${bam}
        modkit pileup --filter-threshold 0.7 --threads ${task.cpus} ${bam} - | awk '\$10 > ${filter}' | bgzip -c > ${bam.simpleName}.bedmethyl.gz
        tabix -p bed ${bam.simpleName}.bedmethyl.gz
        
    else
        echo "No modified reads detected in ${bam}"
        exit 42
    fi

    """
}