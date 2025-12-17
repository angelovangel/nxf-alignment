#!/usr/bin/env nextflow

process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*html'
    
    // some inputs are optional, so better as separated inputs
    input:
        path(runinfo) 
        path(wf_props)
        path(readstats)
        path(ref_stats)
        path(hist) 
        path(bedcov)
        path(bedcov_complement) 
        path(flagstat)
        
    output:
        path "*.html"

    script:
    """
    make-report.py \
    --runinfo $runinfo \
    --wfinfo $wf_props \
    --readstats $readstats \
    --refstats $ref_stats \
    --hist $hist \
    --bedcov $bedcov \
    --bedcov-compl $bedcov_complement \
    --flagstat $flagstat \
    -o nxf-alignment-report.html
    """
}