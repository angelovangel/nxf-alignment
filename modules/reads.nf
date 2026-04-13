#!/usr/bin/env nextflow

process CONVERT_EXCEL {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    input:
        path(excelfile)

    output:
        path("*.csv")

    script:
    """
    convert_excel.R $excelfile
    """
}

// takes in csv, checks for duplicate barcodes etc.
process VALIDATE_SAMPLESHEET {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    //publishDir "$params.outdir", mode: 'copy', pattern: '00-validated-samplesheet.csv'

    input: 
    path(csv)

    output:
    path("samplesheet-validated.csv")
    path "versions.txt", emit: versions

    script:
    """
    validate_samplesheet.R $csv

    cat <<EOF > versions.txt
    ${task.process}: R v\$(R --version | head -n 1 | sed 's/^R version //')
    EOF
    """
}

process MERGE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    errorStrategy {
        if (task.exitStatus == 42) {
            println "MERGE_READS: [WARNING] ${samplename} (${barcode}) defined in the samplesheet but does not exist in the data. Ignoring."
            return 'ignore'
        }
        return 'ignore'
    }
    tag "${samplename} - ${barcode}"

    publishDir "$params.outdir/00-basecall/processed", mode: 'copy', pattern: '*{fastq.gz,fastq,bam}'

    input:
        tuple val(samplename), val(barcode), path(bam_pass)
    
    output: 
        path('*{fastq.gz,fastq,bam}')
    
    script:
    """
    if [ -d "${bam_pass}/${barcode}" ]; then
        samtools cat ${bam_pass}/${barcode}/*.bam > ${samplename}.bam
    else
        exit 42
    fi
    """
}

process READ_STATS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    publishDir "${params.outdir}/00-basecall/readqc", mode: 'copy', pattern: '*readstats.tsv'
    tag "${reads.simpleName}, ${reads.extension} file"

    input:
        path reads

    output:
        path "*readstats.tsv"
        path "versions.txt", emit: versions

    script:
    
    """
    printf "file\\treads\\tbases\\tn_bases\\tmin_len\\tmax_len\\tn50\\tGC_percent\\tQ20_percent\\tmods\\n" > ${reads.simpleName}.readstats.tsv
    faster2 -ts ${reads} | tr -d '\\n' >> ${reads.simpleName}.readstats.tsv
    echo -e "\\t-" >> ${reads.simpleName}.readstats.tsv

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: faster2 v\$(faster2 --version 2>&1 | sed 's/^faster2 //')
    ${task.process}: fasterplot v\$(fasterplot --version 2>&1 | sed 's/^fasterplot //')
    ${task.process}: samtools v\$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """
}

process READ_HIST {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${reads.simpleName}"

    input:
        path reads

    output:
        path "*.hist"

    script:
    """
    faster2 --len ${reads}  | bincount.awk -v type=len  | sort -n > ${reads.simpleName}.len.hist
    faster2 --gc ${reads}   | bincount.awk -v type=gc   | sort -n > ${reads.simpleName}.gc.hist
    fasterplot -q ${reads} | sort -n > ${reads.simpleName}.qual.hist
    """
}
process CONVERT_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    tag "${reads.simpleName}"
    cpus 4

    input:
        path reads

    output:
        path "*.fastq.gz"

    script:
    """
    samtools fastq -@ ${task.cpus} -T '*' ${reads} | gzip > ${reads.simpleName}.fastq.gz
    """
}

process READ_ANI {
    container 'docker.io/staphb/sylph:latest'
    publishDir "${params.outdir}/00-basecall/readqc", mode: 'copy', pattern: '*ani.tsv'
    tag "${reads.simpleName}"

    input:
        tuple path(ref), path(reads)

    output:
        path "*.tsv"

    script:
    """
    sylph query -m 80 ${ref} ${reads} > ${reads.simpleName}.ani.tsv
    """
}
