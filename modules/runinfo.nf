#!/usr/bin/env nextflow

// get run info from bam header
process RUN_INFO {
    
    container 'docker.io/aangeloo/nxf-tgs:latest' 
    
    input:
        // Expects one BAM file as input
        path(bam)

    output:
        // Outputs a CSV file with the extracted info
        path("run_info.csv")

    script:
    """
    
    RG_LINE=\$(samtools view -H ${bam} | grep '^@RG' | head -n 1) 
    
    FLOWCELL_ID=\$(echo \$RG_LINE | sed -n 's/.*PU:\\([^[:space:]]*\\).*/\\1/p')
    BASECALL_MODEL=\$(echo \$RG_LINE | sed -n 's/.*basecall_model=\\([^[:space:]]*\\).*/\\1/p')
    RUN_DATE=\$(echo "\$RG_LINE" | sed -n 's/.*DT:\\([^T]*\\).*/\\1/p')
    RUN_ID=\$(echo "\$RG_LINE" | sed -n 's/.*DS:runid=\\([^[:space:]]*\\).*/\\1/p')

    echo "flowcell_id,basecall_model,run_date,run_id" > run_info.csv
    echo "\$FLOWCELL_ID,\$BASECALL_MODEL,\$RUN_DATE,\$RUN_ID" >> run_info.csv
    """
}
