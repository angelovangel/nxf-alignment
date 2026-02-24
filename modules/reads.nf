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

    script:
    """
    validate_samplesheet.R $csv
    """
}

process MERGE_READS {
    container 'docker.io/aangeloo/nxf-tgs:latest'
    errorStrategy {
        if (task.exitStatus == 42) {
            println "MERGE_READS: [WARNING] ${barcode} defined in the samplesheet but does not exist in the data. Ignoring."
            return 'ignore'
        }
        return 'ignore'
    }
    tag "${barcode} == ${samplename}"

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
    publishDir "${params.outdir}/00-basecall", mode: 'copy', pattern: '*readstats.tsv'
    tag "${reads.simpleName}, ${reads.extension} file"

    input:
        path reads

    output:
        path "*readstats.tsv"
        path "versions.txt", emit: versions

    script:
    
    """
    printf "file\\treads\\tbases\\tn_bases\\tmin_len\\tmax_len\\tn50\\tGC_percent\\tQ20_percent\\tmods\\n" > ${reads.simpleName}.readstats.tsv
    
    mods="-"
    if [[ ${reads.extension} == "bam" ]]; then
        # Detect modifications
        bamtags=\$(samtools view ${reads} | head -n 1000 | grep -o "MM:Z:[^[:space:]]*" | cut -d',' -f1 | cut -d':' -f3 | sort -u)
        mods=\$(echo \$bamtags | parse_mods.py)

        samtools fastq ${reads} | faster2 -ts - | tr -d '\\n' >> ${reads.simpleName}.readstats.tsv
        echo -e "\\t\$mods" >> ${reads.simpleName}.readstats.tsv
    else 
        faster2 -ts ${reads} | tr -d '\\n' >> ${reads.simpleName}.readstats.tsv
        echo -e "\\t-" >> ${reads.simpleName}.readstats.tsv
    fi

    cat <<-END_VERSIONS > versions.txt
    ${task.process}: faster2 v\$(faster2 --version 2>&1 | sed 's/^faster2 //')
    ${task.process}: samtools v\$(samtools --version | head -n 1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
