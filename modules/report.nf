#!/usr/bin/env nextflow

process REPORT {
    container 'docker.io/aangeloo/nxf-tgs:latest'

    publishDir "${params.outdir}", mode: 'copy', pattern: '*html'
    
    input:
       tuple path(hist), path(readstats), path(bedcov), path(bedcov_complement), path(runinfo), path(wf_props), path(ref_stats)
    output:
        path "*.html"

    script:
    """
    make-report.py \
    --hist $hist \
    --readstats $readstats \
    --bedcov $bedcov \
    --bedcov-compl $bedcov_complement \
    --runinfo $runinfo \
    --wfinfo $wf_props \
    --refstats $ref_stats \
    -o nxf-alignment-report.html
    """
}