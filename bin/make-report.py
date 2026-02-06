#!/usr/bin/env python3
"""
Generate an HTML report from coverage histogram files (.hist) and readstats files (.readstats.tsv)
Usage: python bedtools-report.py --hist sample1.hist [sample2.hist ...] --readstats sample1.readstats.tsv [sample2.readstats.tsv ...] -o output.html [--runinfo run_info.csv] [--wfinfo wf_properties.csv]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime
import csv
import json
import gzip

def format_si(num):
    """Format number with SI suffix (K, M, G, T, P)"""
    if num == 0:
        return "0"
    
    suffixes = ['', ' K', ' M', ' G', ' T', ' P']
    suffix_index = 0
    num_float = float(num)
    
    while abs(num_float) >= 1000 and suffix_index < len(suffixes) - 1:
        num_float /= 1000.0
        suffix_index += 1

    return f"{num_float:.1f}{suffixes[suffix_index]}"

def parse_hist_file(filepath):
    """Parse a .hist file and return list of records"""
    data = []
    with open(filepath, 'r') as f:
        # Skip the first line, assuming it is a header
        next(f, None)
        
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 8:
                data.append({
                    'chr': parts[0],
                    'start': int(parts[1]) if parts[1].isdigit() else 0,
                    'end': int(parts[2]) if parts[2].isdigit() else 0,
                    'gene': parts[3],
                    'coverage': int(parts[4]),
                    'count': int(parts[5]),
                    'total': int(parts[6]),
                    'fraction': float(parts[7])
                })
    return data

def parse_readstats_file(filepath):
    """Parse a readstats.tsv file and return a dictionary with stats"""
    with open(filepath, 'r') as f:
        # Read header
        header = next(f).strip().split('\t')
        # Read data line
        data_line = next(f).strip().split('\t')
        
        # Create dictionary from header and data
        stats = {}
        for i, col in enumerate(header):
            if i < len(data_line):
                # Try to convert to appropriate type
                value = data_line[i]
                if col in ['reads', 'bases', 'n_bases', 'min_len', 'max_len', 'n50']:
                    stats[col] = int(value) if value != '-' else 0
                elif col in ['GC_percent', 'Q20_percent']:
                    stats[col] = float(value) if value != '-' else 0.0
                else:
                    stats[col] = value
        
    return stats

def parse_bedcov_file(filepath):
    """Parse a samtools bedcov file and return total length and coverage"""
    total_len = 0
    total_cov = 0
    has_data = False
    
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # chr, start, end, [name], coverage
            # Coverage count is always the last column in default bedcov output
            if len(parts) >= 4:
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    # Last two columns are bases, reads
                    cov = int(parts[-2]) 
                    total_len += (end - start)
                    total_cov += cov
                    has_data = True
                except ValueError:
                    continue
    
    # Return None if file had no valid data
    if not has_data:
        return None
    
    return {'len': total_len, 'cov': total_cov}

def parse_bedcov_per_region(filepath):
    """Parse a samtools bedcov file and return per-region data"""
    regions = []
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                try:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    # If name is present (parts[3]), use it, otherwise use chr:start-end
                    name = parts[3] if len(parts) > 4 else f"{chrom}:{start}-{end}"
                    
                    # Last two columns are always bases, reads
                    cov = int(parts[-2])
                    reads = int(parts[-1])
                    
                    length = end - start
                    mean_cov = cov / length if length > 0 else 0
                    read_len = cov / reads if reads > 0 else 0
                    regions.append({
                        'chr': chrom,
                        'start': start,
                        'end': end,
                        'name': name,
                        'length': length,
                        'bases': cov,
                        'reads': reads,
                        'mean_cov': mean_cov,
                        'read_len': read_len
                    })
                except (ValueError, IndexError):
                    continue
    return regions

def parse_flagstat_file(filepath):
    """Parse a samtools flagstat JSON file"""
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
            # Skip empty JSON or JSON without expected data
            if not data or 'QC-passed reads' not in data:
                return None
            qc_passed = data.get('QC-passed reads', {})
            return {
                'primary_mapped': qc_passed.get('primary mapped', 0),
                'primary_mapped_pct': qc_passed.get('primary mapped %', 0.0)
            }
    except Exception as e:
        print(f"Error parsing flagstat file {filepath}: {e}", file=sys.stderr)
        return None

def parse_bcftools_query(filepath):
    """Parse a bcftools query output file and return variant statistics."""
    stats = {
        'total': 0,
        'snp': 0,
        'indel': 0,
        'pass': 0,
        'high_qual': 0,
        'high_qual_snp': 0,
        'high_qual_indel': 0,
        'ts': 0,
        'tv': 0
    }
    
    transitions = {
        ('A', 'G'), ('G', 'A'),
        ('C', 'T'), ('T', 'C')
    }
    transversions = {
        ('A', 'C'), ('C', 'A'),
        ('A', 'T'), ('T', 'A'),
        ('G', 'C'), ('C', 'G'),
        ('G', 'T'), ('T', 'G')
    }
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 7:
                    continue
                
                # Format: CHROM POS REF ALT TYPE FILTER QUAL [SVTYPE]
                ref = parts[2].upper()
                alt = parts[3].upper()
                var_type = parts[4]
                filter_status = parts[5]
                qual_str = parts[6]

                # Skip non-variants and structural variants (if identified by TYPE or SVTYPE)
                if var_type == 'REF':
                    continue
                
                # If there's an 8th column and it's not empty/., it might be an SV
                if len(parts) >= 8 and parts[7] != '.' and parts[7] != '':
                    continue

                stats['total'] += 1
                
                # Check Quality
                is_high_qual = False
                try:
                    if qual_str != '.' and float(qual_str) >= 30:
                        stats['high_qual'] += 1
                        is_high_qual = True
                except ValueError:
                    pass
                
                if var_type == 'SNP':
                    if is_high_qual: stats['high_qual_snp'] += 1
                elif (var_type == 'INDEL' or var_type == 'SNP,INDEL'):
                    if is_high_qual: stats['high_qual_indel'] += 1

                if filter_status == 'PASS':
                    stats['pass'] += 1
                    
                    if var_type == 'SNP':
                        stats['snp'] += 1
                        # Calculate Ts/Tv for PASS SNPs
                        alts = alt.split(',')
                        for single_alt in alts:
                            pair = (ref, single_alt)
                            if pair in transitions:
                                stats['ts'] += 1
                            elif pair in transversions:
                                stats['tv'] += 1

                    elif (var_type == 'INDEL' or var_type == 'SNP,INDEL'):
                        stats['indel'] += 1
                    
        # Calculate Ts/Tv ratio
        stats['ts_tv_ratio'] = stats['ts'] / stats['tv'] if stats['tv'] > 0 else 0.0
            
        return stats
    except Exception as e:
        print(f"Error parsing query file {filepath}: {e}", file=sys.stderr)
        return None

def parse_sv_query(filepath):
    """Parse an SV query file (from bcftools query) and return statistics."""
    stats = {
        'total': 0,
        'DEL': 0,
        'INS': 0,
        'DUP': 0,
        'INV': 0,
        'BND': 0,
        'other': 0
    }
    
    try:
        with open(filepath, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                
                # Format: CHROM POS REF ALT TYPE FILTER QUAL SVTYPE
                sv_type = parts[7]
                
                if not sv_type or sv_type == '.':
                    # Fallback to ALT if SVTYPE is missing
                    alt = parts[3].strip('<>')
                    if alt in ['DEL', 'INS', 'DUP', 'INV', 'BND']:
                        sv_type = alt
                    else:
                        continue

                stats['total'] += 1
                
                if sv_type in stats:
                    stats[sv_type] += 1
                else:
                    stats['other'] += 1
                    
        return stats
    except Exception as e:
        print(f"Error parsing SV query file {filepath}: {e}", file=sys.stderr)
        return None

def parse_phase_file(filepath):
    """Parse Whatshap phasing stats TSV and return dictionary for 'ALL' chromosome row."""
    try:
        data = {}
        with open(filepath, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                if row['chromosome'] == 'ALL':
                    data = {
                        'phased': int(row.get('phased', 0)),
                        'unphased': int(row.get('unphased', 0)),
                        'blocks': int(row.get('blocks', 0)),
                        'singletons': int(float(row.get('singletons', 0))),
                        'avg_block_bp': float(row.get('bp_per_block_avg', 0)),
                        'phased_fraction': float(row.get('phased_fraction', 0)) * 100
                    }
                    break
        return data if data else None
    except Exception as e:
        print(f"Error parsing phasing stats file {filepath}: {e}", file=sys.stderr)
        return None

def parse_as_file(filepath):
    """Parse adaptive sampling decision file and return accepted/total counts."""
    accepted = 0
    total = 0
    try:
        with open(filepath, 'r') as f:
            # Skip header
            next(f, None)
            for line in f:
                total += 1
                parts = line.strip().split(',')
                if len(parts) >= 2 and parts[1] == 'sequence':
                    accepted += 1
        return accepted, total
    except Exception as e:
        print(f"Error parsing adaptive sampling file {filepath}: {e}", file=sys.stderr)
        return None

def parse_runinfo_csv(filepath):
    """Parse a run info CSV file and return a list of dictionaries (one per row)"""
    run_info = []
    try:
        with open(filepath, 'r', newline='') as f:
            reader = csv.DictReader(f)
            # Read all rows into the list
            for row in reader:
                run_info.append(row)
    except FileNotFoundError:
        print(f"Warning: Run info file not found at {filepath}", file=sys.stderr)
    except Exception as e:
        print(f"Error reading run info CSV: {e}", file=sys.stderr)
        
    return run_info

def calculate_cumulative_coverage(gene_data):
    """Calculate cumulative coverage metrics"""
    total_bases = gene_data[0]['total'] if gene_data else 0
    
    # Calculate bases with coverage >= threshold
    bases_1x = sum(r['count'] for r in gene_data if r['coverage'] >= 1)
    bases_10x = sum(r['count'] for r in gene_data if r['coverage'] >= 10)
    bases_20x = sum(r['count'] for r in gene_data if r['coverage'] >= 20)
    bases_30x = sum(r['count'] for r in gene_data if r['coverage'] >= 30)
    
    return {
        'total': total_bases,
        'cov_1x': bases_1x,
        'cov_10x': bases_10x,
        'cov_20x': bases_20x,
        'cov_30x': bases_30x,
        'pct_1x': (bases_1x / total_bases * 100) if total_bases > 0 else 0,
        'pct_10x': (bases_10x / total_bases * 100) if total_bases > 0 else 0,
        'pct_20x': (bases_20x / total_bases * 100) if total_bases > 0 else 0,
        'pct_30x': (bases_30x / total_bases * 100) if total_bases > 0 else 0,
    }

def get_css():
    """Return the CSS style block by reading from external file"""
    css_path = Path(__file__).parent / 'report' / 'assets' / 'report.css'
    with open(css_path, 'r') as f:
        css_content = f.read()
    return f"""
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/css/select2.min.css" rel="stylesheet" />
    <style>
{css_content}
    </style>
    """

def get_js():
    """Return the Javascript block by reading from external file"""
    js_path = Path(__file__).parent / 'report' / 'assets' / 'report.js'
    with open(js_path, 'r') as f:
        js_content = f.read()
    return f"""
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/js/select2.min.js"></script>
    <script>
{js_content}
    </script>
    """

def render_details_block(title, info_list, add_top_border=False):
    """Render a run info or workflow info details block"""
    if not info_list:
        return ""
    
    border_style = "border-top: 1px solid rgba(255, 255, 255, 0.2);" if add_top_border else ""
    html = f"""
    <details style="padding-top: 5px; margin-top: 5px; {border_style} text-align: left;">
      <summary style="font-size: 0.8em; cursor: pointer; color: white;">{title}</summary>
      <div style="padding-top: 5px; text-align: left; width: 100%;">
        <table style="width: 100%; border-collapse: collapse; margin-top: 5px; color: white; border: none; background: inherit;">
          <tbody>
    """
    
    for i, row in enumerate(info_list):
        if i > 0:
            html += '<tr><td colspan="2" style="border-bottom: 1px solid rgba(255, 255, 255, 0.2); padding: 0;"></td></tr>'
        for key, value in row.items():
            if key is None:
                continue
            display_key = key.replace('_', ' ').title()
            html += f"""
            <tr>
              <td style="padding: 3px 10px 3px 0; font-weight: 500; border: none; background: inherit; font-family: 'Courier New', monospace; color: white; white-space: nowrap; font-size: 0.8em; width: 180px;">{display_key}:</td>
              <td style="padding: 3px 0; border: none; background: inherit; font-family: 'Courier New', monospace; color: white; font-size: 0.8em;">{value}</td>
            </tr>
            """
    
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_stats_cards(readstats_data, samples_data, genes, region_totals, ref_stats, wf_info, as_counts=None):
    """Render the top statistics cards"""
    total_bases = sum(stats.get('bases', 0) for stats in readstats_data.values()) if readstats_data else 0
    total_reads = sum(stats.get('reads', 0) for stats in readstats_data.values()) if readstats_data else 0
    total_bed_size = sum(region_totals.values())
    
    # Adaptive sampling check from wf_info or as_counts
    as_status = "No"
    if as_counts:
        accepted, total = as_counts
        as_status = f"{format_si(accepted)} / {format_si(total)}"
    elif wf_info:
        for row in wf_info:
            if 'Adaptive Sampling' in row:
                as_status = row['Adaptive Sampling']
                break
    
    # Ref stats extraction
    ref_contigs = 0
    ref_bases = 0
    if ref_stats and len(ref_stats) > 0:
        ref_contigs = int(ref_stats[0].get('contigs', 0))
        ref_bases = int(float(ref_stats[0].get('bases', 0))) # in case it's float string

    # Use sorted readstats keys for Sample count definition (consistent with existing)
    sample_count = len(readstats_data) if readstats_data else len(samples_data) # Fallback

    reads_fmt = format_si(total_reads)
    bases_fmt = format_si(total_bases)
    bed_fmt = format_si(total_bed_size)
    ref_bases_fmt = format_si(ref_bases)
    genes_count = len(genes)

    html = f"""
      <div class="stats">
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Number of samples processed in this run.">ⓘ</span> Samples</h3>
          <div class="value">{sample_count}</div>
        </div>
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Total number of sequencing reads across all samples.">ⓘ</span> Total Reads</h3>
          <div class="value">{reads_fmt}</div>
        </div>
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Total number of base pairs sequenced across all samples.">ⓘ</span> Total Bases</h3>
          <div class="value">{bases_fmt}</div>
        </div>
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Total size of the reference genome.">ⓘ</span> Ref Size</h3>
          <div class="value">{ref_bases_fmt}</div>
        </div>
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Number of contigs (chromosomes/scaffolds) in the reference.">ⓘ</span> Ref Contigs</h3>
          <div class="value">{ref_contigs}</div>
        </div>
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Number of regions defined in the target BED file.">ⓘ</span> BED Genes/Regions</h3>
          <div class="value">{genes_count}</div>
        </div>
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Total size of the target regions in the BED file.">ⓘ</span> Total BED size</h3>
          <div class="value">{bed_fmt}</div>
        </div>
        <div class="stat-card">
          <h3><span class="help-icon" data-tooltip="Read statistics (accepted / total reads) from adaptive sampling.">ⓘ</span> Adaptive Sampling</h3>
          <div class="value" style="word-break: break-all;">{as_status}</div>
        </div>
      </div>
    """
    return html

def render_readstats_table(readstats_data):
    """Render the Read Statistics table"""
    if not readstats_data:
        return ""
    
    html = """
      <details class="collapsible-section" open>
        <summary>
          <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
            <h3 style="margin: 0; font-size: 1.1em; color: inherit;">Read Statistics</h3>
            <button class="export-btn" onclick="event.preventDefault(); event.stopPropagation(); exportReadstatsToCSV()">Export CSV</button>
          </div>
        </summary>
        <div class="table-container" style="margin-bottom: 30px; position: relative;">
          <table id="readstatsTable">
            <thead>
              <tr>
                <th class="sample-col sortable" onclick="sortReadstatsTable(0)">Sample</th>
                <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(1)">Reads</th>
                <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(2)">Bases</th>
                
                <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(4)">Min Length</th>
                <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(5)">Max Length</th>
                <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(6)">N50</th>
                <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(7)">GC %</th>
                <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(8)">Q20 %</th>
                <th style="text-align: right;">Mods</th>
              </tr>
            </thead>
            <tbody>
    """
    
    for sample_name in sorted(readstats_data.keys()):
        stats = readstats_data[sample_name]
        html += f"""
            <tr data-sample="{sample_name.lower()}"
                data-reads="{stats.get('reads', 0)}"
                data-bases="{stats.get('bases', 0)}"
                data-nbases="{stats.get('n_bases', 0)}"
                data-minlen="{stats.get('min_len', 0)}"
                data-maxlen="{stats.get('max_len', 0)}"
                data-n50="{stats.get('n50', 0)}"
                data-gc="{stats.get('GC_percent', 0)}"
                data-q20="{stats.get('Q20_percent', 0)}"
                data-mods="{stats.get('mods', '-')}">
              <td class="sample-col">{sample_name}</td>
              <td style="text-align: right;">{stats.get('reads', 0):,}</td>
              <td style="text-align: right;">{stats.get('bases', 0):,}</td>
              
              <td style="text-align: right;">{stats.get('min_len', 0):,}</td>
              <td style="text-align: right;">{stats.get('max_len', 0):,}</td>
              <td style="text-align: right;">{stats.get('n50', 0):,}</td>
              <td style="text-align: right;">{stats.get('GC_percent', 0):.2f}</td>
              <td style="text-align: right;">{stats.get('Q20_percent', 0):.2f}</td>
              <td style="text-align: right;">{stats.get('mods', '-')}</td>
            </tr>
        """
        
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_samtools_table(samtools_data):
    """Render the Samtools Coverage Statistics table"""
    if not samtools_data:
        return ""

    html = """
      <details class="collapsible-section" open>
        <summary>
          <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
            <h3 style="margin: 0; font-size: 1.1em; color: inherit;">Coverage</h3>
            <button class="export-btn" onclick="event.preventDefault(); event.stopPropagation(); exportSamtoolsToCSV()">Export CSV</button>
          </div>
        </summary>
        <div class="table-container" style="margin-bottom: 30px; position: relative;">
          <table id="samtoolsTable">
            <thead>
              <tr>
                <th class="sample-col sortable" onclick="sortSamtoolsTable(0)">Sample</th>
                <th style="text-align: right;" class="sortable" onclick="sortSamtoolsTable(1)">Primary Mapped</th>
                <th style="text-align: right;" class="sortable" onclick="sortSamtoolsTable(2)">Primary Mapped %</th>
                <th style="text-align: right;" class="sortable" onclick="sortSamtoolsTable(3)">Bases on Target</th>
                <th style="text-align: right;" class="sortable" onclick="sortSamtoolsTable(4)">Mean Target Coverage</th>
                <th style="text-align: right;" class="sortable" onclick="sortSamtoolsTable(5)">Bases on Non-target</th>
                <th style="text-align: right;" class="sortable" onclick="sortSamtoolsTable(6)">Mean Non-target Coverage</th>
              </tr>
            </thead>
            <tbody>
    """
    
    for sample_name in sorted(samtools_data.keys()):
        stats = samtools_data[sample_name]
        target = stats.get('target', {'len': 0, 'cov': 0})
        comp = stats.get('non-target', {'len': 0, 'cov': 0})
        flagstat = stats.get('flagstat', {'primary_mapped': 0, 'primary_mapped_pct': 0.0})
        
        target_mean = target['cov'] / target['len'] if target['len'] > 0 else 0
        comp_mean = comp['cov'] / comp['len'] if comp['len'] > 0 else 0

        html += f"""
            <tr data-sample="{sample_name.lower()}"
                data-pmapped="{flagstat['primary_mapped']}"
                data-ppct="{flagstat['primary_mapped_pct']}"
                data-tbases="{target['cov']}"
                data-tcov="{target_mean}"
                data-ntbases="{comp['cov']}"
                data-ntcov="{comp_mean}">
              <td class="sample-col">{sample_name}</td>
              <td style="text-align: right;">{flagstat['primary_mapped']:,}</td>
              <td style="text-align: right;">{flagstat['primary_mapped_pct']:.2f}%</td>
              <td style="text-align: right;">{target['cov']:,}</td>
              <td style="text-align: right;">{target_mean:.2f}</td>
              <td style="text-align: right;">{comp['cov']:,}</td>
              <td style="text-align: right;">{comp_mean:.2f}</td>
            </tr>
        """
        
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_variants_table(variants_data):
    """Render the Variant Statistics table"""
    if not variants_data:
        return ""

    html = """
      <details class="collapsible-section" open>
        <summary>
          <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
            <h3 style="margin: 0; font-size: 1.1em; color: inherit;">Small Variants (SNPs/Indels)</h3>
            <button class="export-btn" onclick="event.preventDefault(); event.stopPropagation(); exportVariantsToCSV()">Export CSV</button>
          </div>
        </summary>
        <div class="table-container" style="margin-bottom: 30px; position: relative;">
          <table id="variantsTable">
            <thead>
              <tr>
                <th class="sample-col sortable" onclick="sortVariantsTable(0)">Sample</th>
                <th style="text-align: right;" class="sortable" onclick="sortVariantsTable(1)">PASS Variants</th>
                <th style="text-align: right;" class="sortable" onclick="sortVariantsTable(2)">SNPs</th>
                <th style="text-align: right;" class="sortable" onclick="sortVariantsTable(3)">Indels</th>
                <th style="text-align: right;" class="sortable" onclick="sortVariantsTable(4)">High Qual (≥Q30)</th>
                <th style="text-align: right;" class="sortable" onclick="sortVariantsTable(5)">High Qual SNPs</th>
                <th style="text-align: right;" class="sortable" onclick="sortVariantsTable(6)">High Qual Indels</th>
                <th style="text-align: right;" class="sortable" onclick="sortVariantsTable(7)">Ts/Tv Ratio</th>
              </tr>
            </thead>
            <tbody>
    """
    
    for sample_name in sorted(variants_data.keys()):
        stats = variants_data[sample_name]
        
        html += f"""
            <tr data-sample="{sample_name.lower()}"
                data-pass="{stats['pass']}"
                data-snp="{stats['snp']}"
                data-indel="{stats['indel']}"
                data-highqual="{stats['high_qual']}"
                data-hqsnp="{stats['high_qual_snp']}"
                data-hqindel="{stats['high_qual_indel']}"
                data-tstv="{stats.get('ts_tv_ratio', 0.0)}">
              <td class="sample-col">{sample_name}</td>
              <td style="text-align: right;">{stats['pass']:,}</td>
              <td style="text-align: right;">{stats['snp']:,}</td>
              <td style="text-align: right;">{stats['indel']:,}</td>
              <td style="text-align: right;">{stats['high_qual']:,}</td>
              <td style="text-align: right;">{stats['high_qual_snp']:,}</td>
              <td style="text-align: right;">{stats['high_qual_indel']:,}</td>
              <td style="text-align: right;">{stats.get('ts_tv_ratio', 0.0):.2f}</td>
            </tr>
        """
        
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_sv_table(sv_data):
    """Render the Structural Variant Statistics table"""
    if not sv_data:
        return ""

    html = """
      <details class="collapsible-section" open>
        <summary>
          <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
            <h3 style="margin: 0; font-size: 1.1em; color: inherit;">Structural Variants (SVs)</h3>
            <button class="export-btn" onclick="event.preventDefault(); event.stopPropagation(); exportSVToCSV()">Export CSV</button>
          </div>
        </summary>
        <div class="table-container" style="margin-bottom: 30px; position: relative;">
          <table id="svTable">
            <thead>
              <tr>
                <th class="sample-col sortable" onclick="sortSVTable(0)">Sample</th>
                <th style="text-align: right;" class="sortable" onclick="sortSVTable(1)">Total SVs</th>
                <th style="text-align: right;" class="sortable" onclick="sortSVTable(2)">Deletions (DEL)</th>
                <th style="text-align: right;" class="sortable" onclick="sortSVTable(3)">Insertions (INS)</th>
                <th style="text-align: right;" class="sortable" onclick="sortSVTable(4)">Duplications (DUP)</th>
                <th style="text-align: right;" class="sortable" onclick="sortSVTable(5)">Inversions (INV)</th>
                <th style="text-align: right;" class="sortable" onclick="sortSVTable(6)">Translocations (BND)</th>
                <th style="text-align: right;" class="sortable" onclick="sortSVTable(7)">Other</th>
              </tr>
            </thead>
            <tbody>
    """
    
    for sample_name in sorted(sv_data.keys()):
        stats = sv_data[sample_name]
        
        html += f"""
            <tr data-sample="{sample_name.lower()}"
                data-total="{stats['total']}"
                data-del="{stats['DEL']}"
                data-ins="{stats['INS']}"
                data-dup="{stats['DUP']}"
                data-inv="{stats['INV']}"
                data-bnd="{stats['BND']}"
                data-other="{stats['other']}">
              <td class="sample-col">{sample_name}</td>
              <td style="text-align: right;">{stats['total']:,}</td>
              <td style="text-align: right;">{stats['DEL']:,}</td>
              <td style="text-align: right;">{stats['INS']:,}</td>
              <td style="text-align: right;">{stats['DUP']:,}</td>
              <td style="text-align: right;">{stats['INV']:,}</td>
              <td style="text-align: right;">{stats['BND']:,}</td>
              <td style="text-align: right;">{stats['other']:,}</td>
            </tr>
        """
        
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_phasing_table(phasing_data):
    """Render the Phasing Statistics table"""
    if not phasing_data:
        return ""

    html = """
      <details class="collapsible-section" open>
        <summary>
          <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
            <h3 style="margin: 0; font-size: 1.1em; color: inherit;">Small Variants (Phasing Statistics)</h3>
            <button class="export-btn" onclick="event.preventDefault(); event.stopPropagation(); exportPhaseToCSV()">Export CSV</button>
          </div>
        </summary>
        <div class="table-container" style="margin-bottom: 30px; position: relative;">
          <table id="phaseTable">
            <thead>
              <tr>
                <th class="sample-col sortable" onclick="sortPhaseTable(0)">Sample</th>
                <th style="text-align: right;" class="sortable" onclick="sortPhaseTable(1)">Phased Variants</th>
                <th style="text-align: right;" class="sortable" onclick="sortPhaseTable(2)">Unphased Variants</th>
                <th style="text-align: right;" class="sortable" onclick="sortPhaseTable(3)">Phased (%)</th>
                <th style="text-align: right;" class="sortable" onclick="sortPhaseTable(4)">Blocks</th>
                <th style="text-align: right;" class="sortable" onclick="sortPhaseTable(5)">Singletons</th>
                <th style="text-align: right;" class="sortable" onclick="sortPhaseTable(6)">Avg Block Size (bp)</th>
              </tr>
            </thead>
            <tbody>
    """
    
    for sample_name in sorted(phasing_data.keys()):
        stats = phasing_data[sample_name]
        
        html += f"""
            <tr data-sample="{sample_name.lower()}"
                data-phased="{stats['phased']}"
                data-unphased="{stats['unphased']}"
                data-pct="{stats['phased_fraction']}"
                data-blocks="{stats['blocks']}"
                data-singletons="{stats['singletons']}"
                data-avgbp="{stats['avg_block_bp']}">
              <td class="sample-col">{sample_name}</td>
              <td style="text-align: right;">{stats['phased']:,}</td>
              <td style="text-align: right;">{stats['unphased']:,}</td>
              <td style="text-align: right;">{stats['phased_fraction']:.2f}%</td>
              <td style="text-align: right;">{stats['blocks']:,}</td>
              <td style="text-align: right;">{stats['singletons']:,}</td>
              <td style="text-align: right;">{stats['avg_block_bp']:.0f}</td>
            </tr>
        """
        
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_coverage_table(samples_data, genes):
    """Render the Coverage Statistics table"""
    # Don't render if there's no data or all samples are empty
    if not samples_data or all(not v for v in samples_data.values()):
        return ""
    
    html = """
      <details class="collapsible-section" open>
        <summary>
          <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
            <h3 style="margin: 0; font-size: 1.1em; color: inherit;">Breadth of Coverage</h3>
            <button class="export-btn" onclick="event.preventDefault(); event.stopPropagation(); exportCoverageToCSV()">Export CSV</button>
          </div>
        </summary>
        <div class="table-container" style="margin-bottom: 30px; position: relative;">
          <table id="dataTable" class="fixed-table">
            <colgroup>
              <col style="width: 180px;">
              <col style="width: 80px;">
              <col style="width: 200px;">
              <col style="width: 120px;">
              <col>
              <col>
              <col>
              <col>
            </colgroup>
            <thead>
              <tr>
                <th rowspan="2" class="sample-col sortable" onclick="sortTable(0)">Sample</th>
                <th rowspan="2" class="chr-col sortable" onclick="sortTable(1)">Chr</th>
                <th rowspan="2" class="gene-col sortable" onclick="sortTable(2)">Gene/Region</th>
                <th rowspan="2" class="size-col sortable" onclick="sortTable(3)">Region size</th>
                <th colspan="4" style="text-align: center; border-bottom: 1px solid #cbd5e1;">Percentage of region with at least X coverage</th>
              </tr>
              <tr>
                <th class="breadth-col sortable" onclick="sortTable(4)">≥1x (%)</th>
                <th class="breadth-col sortable" onclick="sortTable(5)">≥10x (%)</th>
                <th class="breadth-col sortable" onclick="sortTable(6)">≥20x (%)</th>
                <th class="breadth-col sortable" onclick="sortTable(7)">≥30x (%)</th>
              </tr>
            </thead>
            <tbody>
    """
    
    for gene in genes:
        for sample_name, sample_data in samples_data.items():
            gene_data = [r for r in sample_data if r['gene'] == gene]
            if not gene_data:
                continue
                
            location = gene_data[0]
            cov_stats = calculate_cumulative_coverage(gene_data)
            
            html += f"""
            <tr data-gene="{gene.lower()}" data-sample="{sample_name.lower()}" 
                data-chr="{location['chr']}"
                data-total="{cov_stats['total']}" 
                data-cov1="{cov_stats['pct_1x']:.1f}" 
                data-cov10="{cov_stats['pct_10x']:.1f}"
                data-cov20="{cov_stats['pct_20x']:.1f}"
                data-cov30="{cov_stats['pct_30x']:.1f}">
                <td class="sample-col">{sample_name}</td>
                <td class="chr-col">{location['chr']}</td>
                <td class="gene-col"><strong>{gene}</strong></td>
                <td class="size-col">{cov_stats['total']:,}</td>
              <td class="coverage-cell breadth-col">
                <span class="pct-bar" style="width: {cov_stats['pct_1x'] * 0.9}%;">{cov_stats['pct_1x']:.1f}%</span>
              </td>
              <td class="coverage-cell breadth-col">
                <span class="pct-bar" style="width: {cov_stats['pct_10x'] * 0.9}%;">{cov_stats['pct_10x']:.1f}%</span>
              </td>
              <td class="coverage-cell breadth-col">
                <span class="pct-bar" style="width: {cov_stats['pct_20x'] * 0.9}%;">{cov_stats['pct_20x']:.1f}%</span>
              </td>
              <td class="coverage-cell breadth-col">
                <span class="pct-bar" style="width: {cov_stats['pct_30x'] * 0.9}%;">{cov_stats['pct_30x']:.1f}%</span>
              </td>
            </tr>
            """
    
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_bed_coverage_table(bed_coverage_data):
    """Render the Bed Coverage (per-region) table"""
    if not bed_coverage_data:
        return ""

    # Find max mean coverage for scaling bars
    max_mean_cov = 0
    for sample_regions in bed_coverage_data.values():
        for reg in sample_regions:
            if reg['mean_cov'] > max_mean_cov:
                max_mean_cov = reg['mean_cov']
    
    # Ensure we don't divide by zero
    max_mean_cov = max_mean_cov if max_mean_cov > 0 else 1

    html = """
      <details class="collapsible-section" open>
        <summary>
          <div style="display: flex; justify-content: space-between; align-items: center; width: 100%;">
            <h3 style="margin: 0; font-size: 1.1em; color: inherit;">Bed Coverage</h3>
            <button class="export-btn" onclick="event.preventDefault(); event.stopPropagation(); exportBedcovToCSV()">Export CSV</button>
          </div>
        </summary>
        <div class="table-container" style="margin-bottom: 30px; position: relative;">
          <table id="bedcovTable" class="fixed-table">
            <colgroup>
              <col style="width: 180px;">
              <col style="width: 80px;">
              <col style="width: 200px;">
              <col style="width: 120px;">
              <col style="width: 130px;">
              <col style="width: 120px;">
              <col style="width: 120px;">
              <col>
            </colgroup>
            <thead>
              <tr>
                <th class="sample-col sortable" onclick="sortBedcovTable(0)">Sample</th>
                <th class="chr-col sortable" onclick="sortBedcovTable(1)">Chr</th>
                <th class="gene-col sortable" onclick="sortBedcovTable(2)">Gene/Region</th>
                <th class="size-col sortable" onclick="sortBedcovTable(3)">Region Size</th>
                <th class="bases-col sortable" onclick="sortBedcovTable(4)">Bases in Region</th>
                <th class="reads-col sortable" onclick="sortBedcovTable(5)">Reads in Region</th>
                <th class="readlen-col sortable" onclick="sortBedcovTable(6)">Mean Read Len</th>
                <th class="meancov-col sortable" onclick="sortBedcovTable(7)">Mean Coverage</th>
              </tr>
            </thead>
            <tbody>
    """

    for sample_name in sorted(bed_coverage_data.keys()):
        regions = bed_coverage_data[sample_name]
        for reg in regions:
            # Scale bar width to 90% of cell width max
            bar_pct = (reg['mean_cov'] / max_mean_cov) * 90
            html += f"""
                <tr data-sample="{sample_name.lower()}"
                    data-chr="{reg['chr'].lower()}"
                    data-gene="{reg['name'].lower()}"
                    data-start="{reg['start']}"
                    data-end="{reg['end']}"
                    data-length="{reg['length']}"
                    data-bases="{reg['bases']}"
                    data-reads="{reg['reads']}"
                    data-readlen="{reg['read_len']}"
                    data-meancov="{reg['mean_cov']}">
                  <td class="sample-col">{sample_name}</td>
                  <td class="chr-col">{reg['chr']}</td>
                  <td class="gene-col"><strong>{reg['name']}</strong></td>
                  <td class="size-col">{reg['length']:,}</td>
                  <td class="bases-col">{reg['bases']:,}</td>
                  <td class="reads-col">{reg['reads']:,}</td>
                  <td class="readlen-col">{reg['read_len']:.0f}</td>
                  <td class="coverage-cell meancov-col">
                    <span class="pct-bar" style="width: {bar_pct}%;">{reg['mean_cov']:.2f}</span>
                  </td>
                </tr>
            """

    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def generate_html_report(samples_data, readstats_data, run_info, wf_info, ref_stats, samtools_stats, variants_data, sv_data, bed_coverage_data, phasing_data, as_counts, output_file):
    """Generate HTML report from multiple samples"""
    
    # Pre-processing
    all_genes = set()
    for sample_name, sample_data in samples_data.items():
        all_genes.update(r['gene'] for r in sample_data)
    genes = sorted(all_genes)
    
    region_totals = {}
    for gene in genes:
        for sample_name, sample_data in samples_data.items():
            gene_data = [r for r in sample_data if r['gene'] == gene]
            if gene_data:
                region_totals[gene] = gene_data[0]['total']
                break

    # Options for filters
    region_options = "".join([f'<option value="{gene.lower()}">{gene}</option>' for gene in genes])
    
    all_sample_names = set()
    if readstats_data: all_sample_names.update(readstats_data.keys())
    if samples_data: all_sample_names.update(samples_data.keys())
    if samtools_stats: all_sample_names.update(samtools_stats.keys())
    if variants_data: all_sample_names.update(variants_data.keys())
    if bed_coverage_data: all_sample_names.update(bed_coverage_data.keys())
    if phasing_data: all_sample_names.update(phasing_data.keys())
    
    sample_names = sorted(list(all_sample_names))
    sample_options = "".join([f'<option value="{name.lower()}">{name}</option>' for name in sample_names])

    # Render Components
    css_block = get_css()
    js_block = get_js()
    
    run_info_block = render_details_block("Sequencing run details", run_info, add_top_border=True)
    wf_info_block = render_details_block("Workflow details", wf_info, add_top_border=False)
    
    stats_cards = render_stats_cards(readstats_data, samples_data, genes, region_totals, ref_stats, wf_info, as_counts)
    readstats_table = render_readstats_table(readstats_data)
    samtools_table = render_samtools_table(samtools_stats)
    variants_table = render_variants_table(variants_data)
    phasing_table = render_phasing_table(phasing_data)
    sv_table = render_sv_table(sv_data)
    coverage_table = render_coverage_table(samples_data, genes)
    bed_coverage_table = render_bed_coverage_table(bed_coverage_data)
    
    # Assemble HTML
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NXF-ALIGNMENT Report</title>
  {css_block}
</head>
<body>
  <div class="container">
    <div class="header">
      <div class="header-main">
        <h2>NXF-ALIGNMENT Report</h2>
        <a href="https://github.com/angelovangel/nxf-alignment" class="repo-link" target="_blank" title="View Source on GitHub">
          <svg height="24" viewBox="0 0 16 16" version="1.1" width="24" aria-hidden="true" fill="currentColor"><path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path></svg>
        </a>
      </div>
      <p style="text-align: right;">{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
      {run_info_block}
      {wf_info_block} 
    </div>
    <div class="content">
      {stats_cards}
      
      <div class="controls">
        <div class="filter-group">
          <select id="sampleFilter" multiple="multiple" style="width: 100%;">
            {sample_options}
          </select>
        </div>
        <div class="filter-group">
          <select id="regionFilter" multiple="multiple" style="width: 100%;">
            {region_options}
          </select>
        </div>
      </div>
      
      {readstats_table}
      {samtools_table}
      {bed_coverage_table}
      {coverage_table}
      {variants_table}
      {phasing_table}
      {sv_table}
    </div>
  </div>
  <button onclick="scrollToTop()" id="backToTop" title="Go to top">↑</button>
  {js_block}
</body>
</html>
"""

    with open(output_file, 'w') as f:
        f.write(html)
    print(f"HTML report written to {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Generate coverage histogram HTML report')
    parser.add_argument('--hist', nargs='+', default=[], help='One or more .hist.tsv files')
    parser.add_argument('--readstats', nargs='*', default=[], help='One or more .readstats.tsv files')
    parser.add_argument('--runinfo', type=str, help='Optional CSV file with run metadata (e.g., flowcell_id, run_date)', default=None)
    parser.add_argument('--wfinfo', type=str, help='Optional CSV file with workflow properties', default=None)
    parser.add_argument('--refstats', type=str, help='Optional CSV file with reference stats', default=None)
    parser.add_argument('--bedcov', nargs='*', default=[], help='One or more reads.bedcov.tsv files')
    parser.add_argument('--bedcov-compl', nargs='*', default=[], help='One or more reads.bedcov.compl.tsv files')
    parser.add_argument('--flagstat', nargs='*', default=[], help='One or more .flagstat.json files')
    parser.add_argument('--vcf-query', nargs='*', default=[], help='One or more .query files from bcftools query')
    parser.add_argument('--sv-vcf', nargs='*', default=[], help='One or more structural variant VCF files')
    parser.add_argument('--phasestats', nargs='*', default=[], help='One or more phasing stats TSV files')
    parser.add_argument('--asfile', type=str, help='Optional adaptive sampling decision file', default=None)
    parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    
    args = parser.parse_args()
    
    samples_data = {}
    readstats_data = {}
    samtools_stats = {}
    variants_data = {}
    sv_data = {}
    bed_coverage_data = {}
    phasing_data = {}
    run_info = [] 
    wf_info = []
    ref_stats = []
    as_counts = None

    if args.asfile:
        print(f"Parsing adaptive sampling file {args.asfile}...")
        as_counts = parse_as_file(args.asfile)

    if args.runinfo:
        print(f"Loading run info from {args.runinfo}...")
        run_info = parse_runinfo_csv(args.runinfo)
        
    if args.wfinfo:
        print(f"Loading wf info from {args.wfinfo}...")
        wf_info = parse_runinfo_csv(args.wfinfo)

    if args.refstats:
        print(f"Loading ref stats from {args.refstats}...")
        ref_stats = parse_runinfo_csv(args.refstats)

    for hist_file in args.hist:
        path = Path(hist_file)
        sample_name = path.name.replace('.hist.tsv', '')
        # 
        print(f"Processing {hist_file} (Sample: {sample_name})...")
        samples_data[sample_name] = parse_hist_file(hist_file)
        
    for readstats_file in args.readstats:
        path = Path(readstats_file)
        # Assuming filename format is sample.readstats.tsv
        sample_name = path.name.replace('.readstats.tsv', '')
        print(f"Processing stats {readstats_file} (Sample: {sample_name})...")
        readstats_data[sample_name] = parse_readstats_file(readstats_file)

    for query_file in args.vcf_query:
        path = Path(query_file)
        sample_name = path.name.replace('.snp.query', '').replace('.vcf', '').replace('.variants', '')
        print(f"Processing query {query_file} (Sample: {sample_name})...")
        result = parse_bcftools_query(query_file)
        if result is not None:
            variants_data[sample_name] = result

    
    # Actually, sample naming might be simpler: assuming standard nextflow output naming
    # If using replace, we should try to match how sample_name was extracted above.
    
    # Re-iterate carefully.
    for f in args.bedcov:
        # heuristic to remove suffixes
        name = Path(f).name
        # Common suffixes in pipeline
        for suffix in ['.batched.reads.bedcov.tsv', '.reads.bedcov.tsv', '.bedcov.tsv']:
             if name.endswith(suffix):
                 name = name[:-len(suffix)]
                 break
        
        # Summary data
        result = parse_bedcov_file(f)
        if result is not None:
            if name not in samtools_stats: samtools_stats[name] = {}
            samtools_stats[name]['target'] = result
            
        # Per-region data
        regions = parse_bedcov_per_region(f)
        if regions:
            if name not in bed_coverage_data: bed_coverage_data[name] = []
            bed_coverage_data[name].extend(regions)
        
    for f in args.bedcov_compl:
        name = Path(f).name
        for suffix in ['.batched.reads.bedcov.compl.tsv', '.reads.bedcov.compl.tsv', '.bedcov.compl.tsv']:
             if name.endswith(suffix):
                 name = name[:-len(suffix)]
                 break
        result = parse_bedcov_file(f)
        if result is not None:
            if name not in samtools_stats: samtools_stats[name] = {}
            samtools_stats[name]['non-target'] = result

    for f in args.flagstat:
        name = Path(f).name
        # Assuming filename format like sample.flagstat.json
        for suffix in ['.flagstat.json']:
             if name.endswith(suffix):
                 name = name[:-len(suffix)]
                 break
        result = parse_flagstat_file(f)
        if result is not None:
            if name not in samtools_stats: samtools_stats[name] = {}
            samtools_stats[name]['flagstat'] = result

    for vcf_file in args.sv_vcf:
        path = Path(vcf_file)
        sample_name = path.name.replace('.sv.query', '').replace('.query', '').replace('.vcf.gz', '').replace('.vcf', '')
        print(f"Processing SV query {vcf_file} (Sample: {sample_name})...")
        result = parse_sv_query(vcf_file)
        if result is not None:
            sv_data[sample_name] = result
            
    for phase_file in args.phasestats:
        path = Path(phase_file)
        sample_name = path.name.replace('.phase.tsv', '')
        print(f"Processing phasing stats {phase_file} (Sample: {sample_name})...")
        result = parse_phase_file(phase_file)
        if result:
            phasing_data[sample_name] = result

    generate_html_report(samples_data, readstats_data, run_info, wf_info, ref_stats, samtools_stats, variants_data, sv_data, bed_coverage_data, phasing_data, as_counts, args.output)

if __name__ == "__main__":
    main()