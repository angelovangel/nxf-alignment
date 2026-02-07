
process MODKIT {
    container 'docker.io/ontresearch/modkit:latest'
    publishDir "${params.outdir}/04-modifications", mode: 'copy'
    errorStrategy { 
        if (task.exitStatus == 42) {
            println "MODKIT: [WARNING] No modified reads detected for ${task.tag}. Skipping."
            return 'ignore'
        }
        return 'terminate'
    }
    tag "${bam.simpleName}"

    input:
        tuple path(bam), path(bai)

    output:
    path "*.bedmethyl"
    path "*.summary.tsv"
    path "modprobs"

    script:
    """
    modkit modbam summary ${bam} > ${bam.simpleName}.summary.tsv
    nreads=\$(grep "# total_reads_used" ${bam.simpleName}.summary.tsv | awk '{print \$NF}')

    if [ "\$nreads" -gt 0 ]; then
        modkit sample-probs -o modprobs --hist --prefix ${bam.simpleName} ${bam}
        modkit pileup --threads ${task.cpus} ${bam} ${bam.simpleName}.bedmethyl
        
    else
        echo "No modified reads detected in ${bam}"
        exit 42
    fi    
    """
}