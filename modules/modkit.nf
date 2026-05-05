
process MODKIT {
    container 'docker.io/ontresearch/modkit:latest'
    publishDir "${params.outdir}/05-modifications", mode: 'copy'
    errorStrategy { 
        if (task.exitStatus == 42) {
            println "MODKIT: [WARNING] No modified reads detected for ${task.tag}. Skipping."
            return 'ignore'
        }
        return 'ignore'
    }
    tag "${sample}"

    input:
        tuple val(sample), path(bam), path(bai)
        path(fai)

    output:
    path "*.bedmethyl"
    path "*.bw"
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

    modkit modbam summary ${bam} > ${sample}.summary.tsv
    nreads=\$(grep "# total_reads_used" ${sample}.summary.tsv | awk '{print \$NF}')

    if [ "\$nreads" -gt 0 ]; then
        modkit sample-probs -o modprobs --hist --prefix ${sample} ${bam}
        
        # Generate bedmethyl
        modkit pileup --filter-threshold 0.7 --threads ${task.cpus} ${bam} ${sample}.bedmethyl
        
        # Generate BigWig tracks for 5mC (m) and 5hmC (h) if present
        if awk '\$4=="m"' ${sample}.bedmethyl | head -n 1 | grep -q "m"; then
            modkit bedmethyl tobigwig ${sample}.bedmethyl ${sample}.5mC.bw --mod-code m -g ${fai}
        fi
        
        # Try to generate 5hmC if it exists in the data
        if awk '\$4=="h"' ${sample}.bedmethyl | head -n 1 | grep -q "h"; then
            modkit bedmethyl tobigwig ${sample}.bedmethyl ${sample}.5hmC.bw --mod-code h -g ${fai}
        fi
        
    else
        echo "No modified reads detected in ${bam}"
        exit 42
    fi

    """
}