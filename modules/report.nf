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
        path(readhists)
        path(hist) 
        path(bedcov)
        path(bedcov_complement) 
        path(flagstat)
        path(variants)
        path(sv_variants)
        path(phase_stats)
        path(as_file)
        
    output:
        path "*.html"

    script:
    // Check if the file name contains 'empty'
    // If it's empty, we send an empty string; otherwise, we send the full argument
    def hist_arg = hist.name.contains('empty_hist') ? '' : "--hist ${hist}"
    def bedcov_arg = bedcov.name.contains('empty_bedcov') ? '' : "--bedcov ${bedcov}"
    def bedcov_compl_arg = bedcov_complement.name.contains('empty_bedcov_compl') ? '' : "--bedcov-compl ${bedcov_complement}"
    def flagstat_arg = flagstat.name.contains('empty_flagstat') ? '' : "--flagstat ${flagstat}"
    def variants_arg = variants.name.contains('empty_variants') ? '' : "--vcf-query ${variants}"
    def sv_variants_arg = sv_variants.name.contains('empty_sv_variants') ? '' : "--sv-vcf ${sv_variants}"
    def phase_stats_arg = phase_stats.name.contains('empty_phase_stats') ? '' : "--phasestats ${phase_stats}" 
    def asfile_arg = as_file.name.contains('EMPTY') ? '' : "--asfile ${as_file}"
    def readhists_arg = readhists.name.contains('empty_readhists') ? '' : "--readhists ${readhists}"
    """
    make-report.py \
        --runinfo $runinfo \
        --wfinfo $wf_props \
        --readstats $readstats \
        --refstats $ref_stats \
        $hist_arg \
        $bedcov_arg \
        $bedcov_compl_arg \
        $flagstat_arg \
        $variants_arg \
        $sv_variants_arg \
        $phase_stats_arg \
        $phase_stats_arg \
        $asfile_arg \
        $readhists_arg \
        -o nxf-alignment-report.html
    """
}

process VERSIONS {
    publishDir "${params.outdir}/logs", mode: 'copy'
    container 'docker.io/aangeloo/nxf-tgs:latest'
    
    input:
    path 'versions*.txt'
    path summary_file
    
    output:
    path "nxf-alignment-execution-summary.txt"

    script:
    """
    echo -e "\nSoftware Versions" > software_versions.txt
    echo -e "-----------------------------------------" >> software_versions.txt
    
    for f in \$(ls versions*.txt | sort -V); do
        awk -F ': ' '{ if (NF>1) printf "%-25s: %s\\n", \$1, \$2; else print \$0 }' \$f >> software_versions.txt
    done
    
    cat $summary_file software_versions.txt > nxf-alignment-execution-summary.txt
    """
}
