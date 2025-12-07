#!/usr/bin/env python3
"""
Generate an HTML report from coverage histogram files (.hist) and readstats files (.readstats.tsv)
Usage: python bedtools-report.py --hist sample1.hist [sample2.hist ...] --readstats sample1.readstats.tsv [sample2.readstats.tsv ...] -o output.html [--runinfo run_info.csv]
"""

import sys
import argparse
from pathlib import Path
from datetime import datetime
import csv # Import for CSV reading

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


def generate_html_report(samples_data, readstats_data, run_info, output_file):
    """Generate HTML report from multiple samples"""
    
    # Get unique genes
    all_genes = set()
    for sample_name, sample_data in samples_data.items():
        all_genes.update(r['gene'] for r in sample_data)
    genes = sorted(all_genes)
    
    # Pre-generate select options for the Regions filter (used by Select2)
    region_options = "".join([f'<option value="{gene.lower()}">{gene}</option>' for gene in genes])
    
    # Pre-generate select options for the Samples filter (used by Select2)
    sample_names = sorted(samples_data.keys())
    sample_options = "".join([f'<option value="{name.lower()}">{name}</option>' for name in sample_names])
    
    # Calculate total bases per region
    region_totals = {}
    for gene in genes:
        for sample_name, sample_data in samples_data.items():
            gene_data = [r for r in sample_data if r['gene'] == gene]
            if gene_data:
                region_totals[gene] = gene_data[0]['total']
                break
    
    # Calculate global totals for new cards
    total_bases_across_samples = sum(stats.get('bases', 0) for stats in readstats_data.values()) if readstats_data else 0
    total_reads_across_samples = sum(stats.get('reads', 0) for stats in readstats_data.values()) if readstats_data else 0
    total_bed_size = sum(region_totals.values())
    
    # Format numbers using SI suffix
    reads_formatted = format_si(total_reads_across_samples)
    bases_formatted = format_si(total_bases_across_samples)
    bed_formatted = format_si(total_bed_size)

    # Generate run info HTML block with styling matching the .header block and smaller font size
    run_info_html = ""
    if run_info:
        # The <details> block sits inside the .header div, inheriting its dark background.
        run_info_html += """
    <details style="padding-top: 5px; margin-top: 5px; border-top: 1px solid rgba(255, 255, 255, 0.2);">
      <summary style="font-size: 0.8em; cursor: pointer; color: white;">Run Details</summary>
      <div style="padding-top: 5px; text-align: left; max-width: 500px; margin: 0 auto;">
        <table style="width: 100%; border-collapse: collapse; margin-top: 5px; color: white; border: none; background: inherit;">
          <tbody>
    """
        # Assuming only one row of run info for simplicity, but iterating over all
        for i, row in enumerate(run_info):
            if i > 0: # Add separator if there are multiple run info blocks
                run_info_html += '<tr><td colspan="2" style="border-bottom: 1px solid rgba(255, 255, 255, 0.2); padding: 0;"></td></tr>'
            for key, value in row.items():
                # Use clean keys for display
                display_key = key.replace('_', ' ').title()
                run_info_html += f"""
                <tr>
                  <td style="padding: 3px 10px 3px 0; font-weight: 500; border: none; background: inherit; font-family: 'Courier New', monospace; font-weight: 600; color: white; white-space: nowrap; font-size: 0.8em;">{display_key}:</td>
                  <td style="padding: 3px 0; border: none; background: inherit; font-family: 'Courier New', monospace; color: white; font-size: 0.8em;">{value}</td>
                </tr>
"""
        run_info_html += """
            </tbody>
        </table>
    </div>
</details>
"""
    
    # Generate readstats table HTML if data exists
    readstats_table_html = ""
    if readstats_data:
        readstats_table_html = """
      <h3 style="margin-bottom: 15px; font-size: 1.1em; color: #374151;">Read Statistics</h3>
      <div class="table-container" style="margin-bottom: 30px;">
        <table id="readstatsTable">
          <thead>
            <tr>
              <th class="sortable" onclick="sortReadstatsTable(0)">Sample</th>
              <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(1)">Reads</th>
              <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(2)">Bases</th>
              
              <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(4)">Min Length</th>
              <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(5)">Max Length</th>
              <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(6)">N50</th>
              <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(7)">GC %</th>
              <th style="text-align: right;" class="sortable" onclick="sortReadstatsTable(8)">Q20 %</th>
            </tr>
          </thead>
          <tbody>
"""
        
        for sample_name in sorted(readstats_data.keys()):
            stats = readstats_data[sample_name]
            readstats_table_html += f"""
            <tr data-sample="{sample_name.lower()}"
                data-reads="{stats.get('reads', 0)}"
                data-bases="{stats.get('bases', 0)}"
                data-nbases="{stats.get('n_bases', 0)}"
                data-minlen="{stats.get('min_len', 0)}"
                data-maxlen="{stats.get('max_len', 0)}"
                data-n50="{stats.get('n50', 0)}"
                data-gc="{stats.get('GC_percent', 0)}"
                data-q20="{stats.get('Q20_percent', 0)}">
              <td class="sample-col">{sample_name}</td>
              <td style="text-align: right;">{stats.get('reads', 0):,}</td>
              <td style="text-align: right;">{stats.get('bases', 0):,}</td>
              
              <td style="text-align: right;">{stats.get('min_len', 0):,}</td>
              <td style="text-align: right;">{stats.get('max_len', 0):,}</td>
              <td style="text-align: right;">{stats.get('n50', 0):,}</td>
              <td style="text-align: right;">{stats.get('GC_percent', 0):.2f}</td>
              <td style="text-align: right;">{stats.get('Q20_percent', 0):.2f}</td>
            </tr>
"""
        
        readstats_table_html += """
          </tbody>
        </table>
      </div>
      
      <h3 style="margin-bottom: 15px; font-size: 1.1em; color: #374151;">Coverage Statistics</h3>
"""
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>NXF-ALIGNMENT Report</title>
  
  <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
  
  <link href="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/css/select2.min.css" rel="stylesheet" />
  
  <style>
    * {{ margin: 0; padding: 0; box-sizing: border-box; }}
    body {{ 
      /* Use Inter as the main UI font (loaded via CDN) */
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
      background: #f5f5f5;
      min-height: 100vh;
      padding: 10px;
    }}
    .container {{
      max-width: 1800px;
      margin: 0 auto;
      background: white;
      border-radius: 8px;
      box-shadow: 0 2px 8px rgba(0,0,0,0.1);
      overflow: hidden;
    }}
    .header {{
      /* Header color: Dark Gray/Blue */
      background: #374151; 
      color: white;
      padding: 8px 10px; 
      text-align: center;
    }}
    .header h2 {{ 
      font-size: 1.2em; 
      margin-bottom: 2px; 
    }}
    .header p {{ 
      opacity: 0.9; 
      font-size: 0.8em; 
      margin-bottom: 0; /* Ensures the run info block starts right after */
    }}
    .content {{ padding: 10px; }}
    .stats {{
      display: grid;
      /* Adjusted grid template to fit more cards */
      grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); 
      gap: 15px;
      margin-top: 20px; /* Adjusted margin to accommodate run info */
      margin-bottom: 20px;
    }}
    .stat-card {{
      /* CHANGED: Use header color */
      background: #374151; 
      color: white;
      padding: 10px; /* Reduced padding for smaller size */
      border-radius: 8px;
      text-align: center;
    }}
    .stat-card h3 {{ font-size: 0.8em; opacity: 0.9; margin-bottom: 3px; }}
    .stat-card .value {{ 
      font-size: 1.2em; /* Reduced font size */
      font-weight: bold; 
    }}
    .controls {{
      display: flex;
      gap: 20px;
      align-items: flex-end;
      margin-bottom: 20px;
    }}
    .filter-group {{
      flex-grow: 1;
    }}
    .filter-label {{
      display: block;
      margin-bottom: 5px;
      font-weight: 600;
      font-size: 0.9em;
    }}
    .filter-box {{
      padding: 10px;
      width: 100%;
      border: 2px solid #e2e8f0;
      border-radius: 6px;
      font-size: 16px;
    }}
    .filter-box:focus {{
      outline: none;
      border-color: #60a5fa;
    }}
    /* --- Select2 Custom Styling --- */
    .filter-group .select2-container {{
      width: 100% !important; 
    }}
    .select2-container .select2-selection--multiple {{
        min-height: 40px; 
        border: 2px solid #e2e8f0 !important;
        border-radius: 6px;
        padding: 5px 10px;
    }}
    .select2-container--default.select2-container--focus .select2-selection--multiple {{
        border-color: #60a5fa !important;
        outline: 0;
    }}
    .select2-container--default .select2-selection--multiple .select2-selection__choice {{
        background-color: #e0f2fe;
        border: 1px solid #90cdf4;
        color: #1e40af;
        padding: 3px 5px;
        border-radius: 4px;
    }}
    /* --- End Select2 Custom Styling --- */
    
    .table-container {{
      overflow-x: auto;
      border-radius: 8px;
      border: 1px solid #e2e8f0;
      max-height: 600px;
      overflow-y: auto;
      position: relative;
    }}
    table {{
      width: 100%;
      border-collapse: collapse;
      background: white;
    }}
    thead {{
      position: sticky;
      top: 0;
      z-index: 100;
      background: #f8fafc;
    }}
    th, td {{
      padding: 10px;
      text-align: left;
      border-bottom: 1px solid #e2e8f0;
      font-size: 0.9em;
      /* Monospace font for data */
      font-family: 'Courier New', monospace; 
    }}
    th {{
      background: #f8fafc;
      font-weight: 600;
      /* Sans-serif for headers */
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
    }}
    /*tr:hover {{ background: #f8fafc; }}*/
    .coverage-cell {{
      text-align: right;
    }}
    .pct-bar {{
      /* Match card and title background color with reduced opacity */
      display: inline-block;
      height: 16px;
      background: rgba(55, 65, 81, 0.5); 
      border-radius: 2px;
      min-width: 2px;
      vertical-align: middle;
      margin-left: 5px;
    }}
    .sample-col {{
      font-weight: 600;
      background: #f8fafc;
    }}
    .export-btn {{
      /* CHANGED: Use header color */
      background: #374151; 
      color: white;
      border: none;
      padding: 10px 20px;
      border-radius: 6px;
      cursor: pointer;
      font-size: 14px;
      font-weight: 600;
      white-space: nowrap; 
    }}
    .export-btn:hover {{
      /* CHANGED: Darker shade of header color for hover */
      background: #1f2937;
    }}
    th.sortable {{
      cursor: pointer;
      user-select: none;
      position: relative;
      padding-right: 20px;
    }}
    th.sortable:hover {{
      background: #e2e8f0;
    }}
    th.sortable::after {{
      content: '⇅';
      position: absolute;
      right: 5px;
      opacity: 0.3;
      font-size: 0.8em;
    }}
    th.sortable.asc::after {{
      content: '▲';
      opacity: 1;
    }}
    th.sortable.desc::after {{
      content: '▼';
      opacity: 1;
    }}
  </style>
</head>
<body>
  <div class="container">
    <div class="header">
      <h2>NXF-ALIGNMENT Report</h2>
      <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
      {run_info_html} </div>
    <div class="content">
      <div class="stats">
        <div class="stat-card">
          <h3>Total Reads</h3>
          <div class="value">{reads_formatted}</div>
        </div>
        <div class="stat-card">
          <h3>Total Bases</h3>
          <div class="value">{bases_formatted}</div>
        </div>
        <div class="stat-card">
          <h3>Total Genes/Regions</h3>
          <div class="value">{len(genes)}</div>
        </div>
        <div class="stat-card">
          <h3>Samples</h3>
          <div class="value">{len(samples_data)}</div>
        </div>
        <div class="stat-card">
          <h3>Total BED size</h3>
          <div class="value">{bed_formatted}</div>
        </div>
      </div>
      
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
        
        <button class="export-btn" onclick="exportReadstatsToCSV()">Export Read Summary</button>
        <button class="export-btn" onclick="exportCoverageToCSV()">Export Coverage</button>
      </div>
      
      {readstats_table_html}
      
      <div class="table-container">
        <table id="dataTable">
          <thead>
            <tr>
              <th rowspan="2" class="sortable" onclick="sortTable(0)">Sample</th>
              <th rowspan="2" class="sortable" onclick="sortTable(1)">Chr</th>
              <th rowspan="2" class="sortable" onclick="sortTable(2)">Gene/Region</th>
              <th rowspan="2">Location</th>
              <th rowspan="2" class="sortable" onclick="sortTable(4)">Region size</th>
              <th colspan="4" style="text-align: center; border-bottom: 1px solid #cbd5e1;">Percentage of region with at least X coverage</th>
            </tr>
            <tr>
              <th class="sortable" onclick="sortTable(5)" style="text-align: right;">≥1x (%)</th>
              <th class="sortable" onclick="sortTable(6)" style="text-align: right;">≥10x (%)</th>
              <th class="sortable" onclick="sortTable(7)" style="text-align: right;">≥20x (%)</th>
              <th class="sortable" onclick="sortTable(8)" style="text-align: right;">≥30x (%)</th>
            </tr>
          </thead>
          <tbody>
"""
    
    # Generate table rows
    for gene in genes:
        for sample_name, sample_data in samples_data.items():
            gene_data = [r for r in sample_data if r['gene'] == gene]
            
            if not gene_data:
                continue
                
            location = gene_data[0]
            cov_stats = calculate_cumulative_coverage(gene_data)
            
            # Row data attributes are used for sorting, filtering, and clean CSV export
            html += f"""
            <tr data-gene="{gene.lower()}" data-sample="{sample_name.lower()}" 
                data-chr="{location['chr']}"
                data-location="{location['start']:,}-{location['end']:,}"
                data-total="{cov_stats['total']}" 
                data-cov1="{cov_stats['pct_1x']:.1f}" 
                data-cov10="{cov_stats['pct_10x']:.1f}"
                data-cov20="{cov_stats['pct_20x']:.1f}"
                data-cov30="{cov_stats['pct_30x']:.1f}">
              <td class="sample-col">{sample_name}</td>
              <td>{location['chr']}</td>
              <td><strong>{gene}</strong></td>
              <td>{location['start']:,}-{location['end']:,}</td>
              <td style="text-align: right;">{cov_stats['total']:,}</td>
              <td class="coverage-cell">
                {cov_stats['pct_1x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_1x']*1.5}px;"></span>
              </td>
              <td class="coverage-cell">
                {cov_stats['pct_10x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_10x']*1.5}px;"></span>
              </td>
              <td class="coverage-cell">
                {cov_stats['pct_20x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_20x']*1.5}px;"></span>
              </td>
              <td class="coverage-cell">
                {cov_stats['pct_30x']:.1f}%
                <span class="pct-bar" style="width: {cov_stats['pct_30x']*1.5}px;"></span>
              </td>
            </tr>
"""
    
    html += """
          </tbody>
        </table>
      </div>
    </div>
  </div>
  
  <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/js/select2.min.js"></script>
  
  <script>
    let sortDirection = {};
    let readstatsSortDirection = {};
    
    // Utility function to trigger CSV download
    function downloadCSV(csvContent, filename) {
        const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
        const link = document.createElement('a');
        const url = URL.createObjectURL(blob);
        link.setAttribute('href', url);
        link.setAttribute('download', filename);
        link.style.visibility = 'hidden';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }

    // NEW FUNCTION for Read Stats Export
    function exportReadstatsToCSV() {
      const table = document.getElementById('readstatsTable');
      if (!table) {
        alert('Read Statistics table not found.');
        return;
      }

      let csv = [];
      
      // Headers (mirroring the table columns)
      const headers = [
        'Sample', 
        'Reads', 
        'Bases', 
        'Min_Length', 
        'Max_Length', 
        'N50', 
        'GC_Percent', 
        'Q20_Percent'
      ];
      csv.push(headers.join(','));
      
      // Data rows (only visible ones)
      const dataRows = table.querySelectorAll('tbody tr');
      dataRows.forEach(row => {
        if (row.style.display !== 'none') {
          // Use data attributes for clean data extraction
          const cols = [
            row.getAttribute('data-sample'),
            row.getAttribute('data-reads'),
            row.getAttribute('data-bases'),
            row.getAttribute('data-minlen'),
            row.getAttribute('data-maxlen'),
            row.getAttribute('data-n50'),
            row.getAttribute('data-gc'),
            row.getAttribute('data-q20'),
          ];
          
          // Escape quotes and wrap in quotes if contains comma
          const safeCols = cols.map(text => {
            text = text.trim();
            // Remove commas from numbers
            text = text.replace(/,/g, ''); 
            
            if (text.includes('"')) {
              text = '"' + text.replace(/"/g, '""') + '"';
            }
            return text;
          });
          
          csv.push(safeCols.join(','));
        }
      });
      
      const csvContent = csv.join('\\n');
      downloadCSV(csvContent, 'readstats_report_export.csv');
    }
    
    // RENAMED exportToCSV to exportCoverageToCSV
    function exportCoverageToCSV() {
      const table = document.getElementById('dataTable');
      let csv = [];
      
      // Headers
      const headers = [
        'Sample', 
        'Chr', 
        'Gene/Region', 
        'Location', 
        'Region_Size', 
        'Pct_Ge_1x', 
        'Pct_Ge_10x', 
        'Pct_Ge_20x', 
        'Pct_Ge_30x'
      ];
      csv.push(headers.join(','));
      
      // Data rows (only visible ones)
      const dataRows = table.querySelectorAll('tbody tr');
      dataRows.forEach(row => {
        if (row.style.display !== 'none') {
          // Use data attributes for clean data extraction
          const cols = [
            row.getAttribute('data-sample'),
            row.getAttribute('data-chr'),
            row.getAttribute('data-gene'),
            // Keep original location format (e.g., "1,000-2,000") but remove commas for clean CSV
            row.getAttribute('data-location').replace(/,/g, ''), 
            row.getAttribute('data-total'),
            row.getAttribute('data-cov1'),
            row.getAttribute('data-cov10'),
            row.getAttribute('data-cov20'),
            row.getAttribute('data-cov30'),
          ];
          
          // Escape quotes and wrap in quotes if contains comma
          const safeCols = cols.map(text => {
            text = text.trim();
            if (text.includes(',') || text.includes('"')) {
              text = '"' + text.replace(/"/g, '""') + '"';
            }
            return text;
          });
          
          csv.push(safeCols.join(','));
        }
      });
      
      const csvContent = csv.join('\\n');
      downloadCSV(csvContent, 'coverage_report_export.csv');
    }

    
    // --- Select2 Initialization and Filter Binding ---
    $(document).ready(function() {
      // Initialize Select2 on the sample filter
      $('#sampleFilter').select2({
        placeholder: "Filter samples...",
        allowClear: true,
        closeOnSelect: true 
      }).on('change', function() {
          // Trigger filtering for both tables on sample change
          filterTable();
          filterReadstatsTable(); 
      }); 
      
      // Initialize Select2 on the region filter
      $('#regionFilter').select2({
        placeholder: "Filter regions...",
        allowClear: true,
        closeOnSelect: true 
      }).on('change', filterTable); // Only Coverage table uses region filter
      
      // Apply initial filter
      filterTable();
      filterReadstatsTable(); 
      
      // Default Sort: Sort both tables by Sample Name (Column 0) on load.
      // Set initial direction to 'desc' so the toggle makes it 'asc'
      sortDirection[0] = 'desc'; 
      readstatsSortDirection[0] = 'desc';
      
      // Call sort to apply default 'asc' sort
      sortTable(0);
      sortReadstatsTable(0);
    });
    // -------------------------------------------------
    
    function sortReadstatsTable(columnIndex) {
      const table = document.getElementById('readstatsTable');
      if (!table) return;
      
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.querySelectorAll('tr'));
      
      // Toggle sort direction
      readstatsSortDirection[columnIndex] = readstatsSortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
      const isAsc = readstatsSortDirection[columnIndex] === 'asc';
      
      // Update header styling
      const headers = table.querySelectorAll('th.sortable');
      headers.forEach(h => {
        h.classList.remove('asc', 'desc');
      });
      // Check if the current column exists and is a header before adding classes
      if (headers[columnIndex]) {
        headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');
      }
      
      rows.sort((a, b) => {
        let aVal, bVal;
        
        // Corrected index map for readstats columns to data attributes
        const dataAttrMap = {
            0: 'sample',
            1: 'reads',
            2: 'bases',
            // Note: Indices 3 are missing in HTML (n_bases), so indices 4,5,6,7,8 map to min_len onwards
            4: 'minlen', // Column 4 in table
            5: 'maxlen', // Column 5
            6: 'n50',    // Column 6
            7: 'gc',     // Column 7
            8: 'q20'     // Column 8
        };
        
        const dataKey = dataAttrMap[columnIndex];

        if (!dataKey) {
            // Treat unmapped columns as non-sortable or use default string sort for safety
             aVal = a.cells[columnIndex].textContent.toLowerCase().trim();
             bVal = b.cells[columnIndex].textContent.toLowerCase().trim();
             return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        }

        aVal = a.getAttribute('data-' + dataKey);
        bVal = b.getAttribute('data-' + dataKey);
        
        if (columnIndex === 0) { // Sample name (string)
          return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        } else {
          aVal = parseFloat(aVal);
          bVal = parseFloat(bVal);
          return isAsc ? aVal - bVal : bVal - aVal;
        }
      });
      
      rows.forEach(row => tbody.appendChild(row));
    }

    function filterReadstatsTable() {
        const selectedSamples = $('#sampleFilter').val() || [];
        const table = document.getElementById('readstatsTable');
        if (!table) return;

        const rows = table.getElementsByTagName('tbody')[0].getElementsByTagName('tr');

        for (let i = 0; i < rows.length; i++) {
            const row = rows[i];
            const sample = row.getAttribute('data-sample');
            
            // If no samples are selected, show all. Otherwise, only show selected.
            let sampleMatch = true;
            if (selectedSamples.length > 0) {
                sampleMatch = selectedSamples.includes(sample);
            }
            
            row.style.display = sampleMatch ? '' : 'none';
        }
    }
    
    function sortTable(columnIndex) {
      const table = document.getElementById('dataTable');
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.querySelectorAll('tr'));
      
      // Toggle sort direction
      sortDirection[columnIndex] = sortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
      const isAsc = sortDirection[columnIndex] === 'asc';
      
      // Update header styling
      const headers = table.querySelectorAll('th.sortable');
      headers.forEach(h => {
        h.classList.remove('asc', 'desc');
      });
      // Check if the current column exists and is a header before adding classes
      if (headers[columnIndex]) {
        headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');
      }
      
      rows.sort((a, b) => {
        let aVal, bVal;
        
        if (columnIndex === 0) { // Sample
          aVal = a.getAttribute('data-sample');
          bVal = b.getAttribute('data-sample');
        } else if (columnIndex === 1) { // Chr
          aVal = a.getAttribute('data-chr');
          bVal = b.getAttribute('data-chr');
        } else if (columnIndex === 2) { // Gene
          aVal = a.getAttribute('data-gene');
          bVal = b.getAttribute('data-gene');
        } else if (columnIndex === 3) { // Location (Not sortable - fallback logic not needed since it's not marked sortable in HTML)
          // Fallback to sorting by Sample then Gene for stable order
          if (a.getAttribute('data-sample') !== b.getAttribute('data-sample')) {
            return isAsc ? a.getAttribute('data-sample').localeCompare(b.getAttribute('data-sample')) : b.getAttribute('data-sample').localeCompare(a.getAttribute('data-sample'));
          }
          aVal = a.getAttribute('data-gene');
          bVal = b.getAttribute('data-gene');
        } else if (columnIndex === 4) { // Region Size (data-total)
          aVal = parseFloat(a.getAttribute('data-total'));
          bVal = parseFloat(b.getAttribute('data-total'));
        } else if (columnIndex === 5) { // ≥1x (data-cov1)
          aVal = parseFloat(a.getAttribute('data-cov1'));
          bVal = parseFloat(b.getAttribute('data-cov1'));
        } else if (columnIndex === 6) { // ≥10x (data-cov10)
          aVal = parseFloat(a.getAttribute('data-cov10'));
          bVal = parseFloat(b.getAttribute('data-cov10'));
        } else if (columnIndex === 7) { // ≥20x (data-cov20)
          aVal = parseFloat(a.getAttribute('data-cov20'));
          bVal = parseFloat(b.getAttribute('data-cov20'));
        } else if (columnIndex === 8) { // ≥30x (data-cov30)
          aVal = parseFloat(a.getAttribute('data-cov30'));
          bVal = parseFloat(b.getAttribute('data-cov30'));
        }
        
        if (typeof aVal === 'string') {
          return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
        } else {
          // Handle cases where aVal or bVal might be NaN/null for a safer comparison
          if (isNaN(aVal) || isNaN(bVal)) {
            return 0;
          }
          return isAsc ? aVal - bVal : bVal - aVal;
        }
      });
      
      rows.forEach(row => tbody.appendChild(row));
    }
    
    function filterTable() {
      // Get selected samples from Select2 (returns an array of values, or null if none selected)
      const selectedSamples = $('#sampleFilter').val();
      const sampleFilters = selectedSamples || [];
      
      // Get selected regions from Select2
      const selectedRegions = $('#regionFilter').val(); 
      const regionFilters = selectedRegions || [];
      
      const table = document.getElementById('dataTable');
      const rows = table.getElementsByTagName('tbody')[0].getElementsByTagName('tr');
      
      for (let i = 0; i < rows.length; i++) {
        const row = rows[i];
        const sample = row.getAttribute('data-sample');
        const gene = row.getAttribute('data-gene');
        
        // 1. Sample Match: If filters exist, the row sample must be in the selected list.
        let sampleMatch = true;
        if (sampleFilters.length > 0) {
            sampleMatch = sampleFilters.includes(sample);
        }
        
        // 2. Region Match: If filters exist, the row gene must be in the selected list.
        let regionMatch = true;
        if (regionFilters.length > 0) {
          regionMatch = regionFilters.includes(gene);
        }
        
        if (sampleMatch && regionMatch) {
          row.style.display = '';
        } else {
          row.style.display = 'none';
        }
      }
    }
  </script>
</body>
</html>
"""
    
    with open(output_file, 'w') as f:
        f.write(html)
    
    print(f"Report generated: {output_file}")


def main():
    parser = argparse.ArgumentParser(description='Generate coverage histogram HTML report')
    parser.add_argument('--hist', nargs='+', required=True, help='One or more .hist files')
    parser.add_argument('--readstats', nargs='*', default=[], help='One or more .readstats.tsv files')
    # NEW ARGUMENT
    parser.add_argument('--runinfo', type=str, help='Optional CSV file with run metadata (e.g., flowcell_id, run_date)', default=None)
    parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    
    args = parser.parse_args()
    
    samples_data = {}
    readstats_data = {}
    run_info = [] # Initialize as empty list

    # Load run info if provided
    if args.runinfo:
        print(f"Loading run info from {args.runinfo}...")
        run_info = parse_runinfo_csv(args.runinfo)

    for sample_file in args.hist:
        sample_name = Path(sample_file).stem.replace('.hist', '')
        print(f"Loading {sample_file}...")
        samples_data[sample_name] = parse_hist_file(sample_file)
    
    for readstats_file in args.readstats:
        sample_name = Path(readstats_file).stem.replace('.readstats', '')
        print(f"Loading {readstats_file}...")
        readstats_data[sample_name] = parse_readstats_file(readstats_file)
    
    print("Generating HTML report...")
    # Pass the run_info list to the generate function
    generate_html_report(samples_data, readstats_data, run_info, args.output)
    
    print(f"\nDone! Open {args.output} in your browser to view the report.")
    print(f"Total samples: {len(samples_data)}")
    print(f"Sample names: {', '.join(samples_data.keys())}")
    if readstats_data:
        print(f"Readstats available for: {', '.join(readstats_data.keys())}")

if __name__ == "__main__":
    main()