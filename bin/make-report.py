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
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            # chr, start, end, [name], coverage
            # Coverage count is always the last column in default bedcov output
            if len(parts) >= 4:
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    cov = int(parts[-1]) 
                    total_len += (end - start)
                    total_cov += cov
                except ValueError:
                    continue
    return {'len': total_len, 'cov': total_cov}

def parse_flagstat_file(filepath):
    """Parse a samtools flagstat JSON file"""
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
            qc_passed = data.get('QC-passed reads', {})
            return {
                'primary_mapped': qc_passed.get('primary mapped', 0),
                'primary_mapped_pct': qc_passed.get('primary mapped %', 0.0)
            }
    except Exception as e:
        print(f"Error parsing flagstat file {filepath}: {e}", file=sys.stderr)
        return {'primary_mapped': 0, 'primary_mapped_pct': 0.0}

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
    """Return the CSS style block"""
    return """
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <link href="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/css/select2.min.css" rel="stylesheet" />
    <style>
      * { margin: 0; padding: 0; box-sizing: border-box; }
      body { 
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
        background: #f5f5f5;
        min-height: 100vh;
        padding: 10px;
      }
      .container {
        max-width: 1800px;
        margin: 0 auto;
        background: white;
        border-radius: 8px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.1);
        overflow: hidden;
      }
      .header {
        background: #374151; 
        color: white;
        padding: 8px 10px; 
        text-align: center;
      }
      .header h2 { 
        font-size: 1.2em; 
        margin-bottom: 2px; 
      }
      .header p { 
        opacity: 0.9; 
        font-size: 0.8em; 
        margin-bottom: 0; 
      }
      .content { padding: 10px; }
      .stats {
        display: grid;
        grid-template-columns: repeat(auto-fit, minmax(140px, 1fr)); 
        gap: 15px;
        margin-top: 0px; 
        margin-bottom: 10px;
      }
      .stat-card {
        background: #374151; 
        color: white;
        padding: 10px; 
        border-radius: 8px;
        text-align: center;
      }
      .stat-card h3 { font-size: 0.8em; opacity: 0.9; margin-bottom: 3px; }
      .stat-card .value { 
        font-size: 1.2em; 
        font-weight: bold; 
      }
      .controls {
        display: flex;
        gap: 20px;
        align-items: flex-end;
        margin-bottom: 20px;
      }
      .filter-group { flex-grow: 1; }
      .filter-label {
        display: block;
        margin-bottom: 5px;
        font-weight: 400;
        font-size: 0.9em;
      }
      .filter-box {
        padding: 10px;
        width: 100%;
        border: 2px solid #374151;
        background-color: #374151;
        border-radius: 6px;
        font-size: 16px;
      }
      .filter-box:focus { outline: none; border-color: #374151; }
      .filter-group .select2-container { width: 'resolve' !important; }
      .table-container {
        overflow-x: auto;
        border-radius: 8px;
        border: 1px solid #e2e8f0;
        max-height: 600px;
        overflow-y: auto;
        position: relative;
      }
      table {
        width: 100%;
        border-collapse: collapse;
        background: white;
      }
      thead {
        position: sticky;
        top: 0;
        z-index: 100;
        background: #f8fafc;
      }
      th, td {
        padding: 10px;
        text-align: left;
        border-bottom: 1px solid #e2e8f0;
        font-size: 0.9em;
        font-family: 'Courier New', monospace; 
      }
      th {
        background: #f8fafc;
        font-weight: 600;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', Arial, sans-serif;
      }
      .coverage-cell { text-align: left; }
      .pct-bar {
        display: block;
        height: 16px;
        background: rgba(55, 65, 81, 0.3); 
        border-radius: 2px;
        min-width: 2px;
        vertical-align: middle;
        margin-left: 20px;
      }
      .pct-bar:hover { background: rgba(55, 65, 81, 0.5); }
      .sample-col { 
        font-weight: 600; 
        background: #f8fafc;
        width: 100px;
        min-width: 100px;
        max-width: 200px;
        white-space: nowrap;
        overflow: hidden;
        text-overflow: ellipsis;
      }
      .export-btn {
        background: #ffffff; 
        color: #6b7280;
        border: 1px solid #e5e7eb;
        padding: 5px 10px;
        border-radius: 4px;
        cursor: pointer;
        font-size: 0.85em;
        transition: all 0.2s;
        font-weight: 600;
      }
      .export-btn:hover { 
        background: #f9fafb; 
        border-color: #d1d5db; 
        color: #374151; 
      }
      th.sortable {
        cursor: pointer;
        user-select: none;
        position: relative;
        padding-right: 20px;
      }
      th.sortable:hover { background: #e2e8f0; }
      th.sortable::after {
        content: '⇅';
        position: absolute;
        right: 5px;
        opacity: 0.3;
        font-size: 0.8em;
      }
      th.sortable.asc::after { content: '▲'; opacity: 1; }
      th.sortable.desc::after { content: '▼'; opacity: 1; }
    </style>
    """

def get_js():
    """Return the Javascript block"""
    return """
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/select2/4.0.13/js/select2.min.js"></script>
    <script>
        let sortDirection = {};
        let readstatsSortDirection = {};
        let samtoolsSortDirection = {};
        
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

        function exportReadstatsToCSV() {
          const table = document.getElementById('readstatsTable');
          if (!table) { alert('Read Statistics table not found.'); return; }
          let csv = [];
          const headers = ['Sample', 'Reads', 'Bases', 'Min_Length', 'Max_Length', 'N50', 'GC_Percent', 'Q20_Percent'];
          csv.push(headers.join(','));
          const dataRows = table.querySelectorAll('tbody tr');
          dataRows.forEach(row => {
            if (row.style.display !== 'none') {
              const cols = [
                row.getAttribute('data-sample'), row.getAttribute('data-reads'), row.getAttribute('data-bases'),
                row.getAttribute('data-minlen'), row.getAttribute('data-maxlen'), row.getAttribute('data-n50'),
                row.getAttribute('data-gc'), row.getAttribute('data-q20')
              ];
              const safeCols = cols.map(text => {
                text = text.trim().replace(/,/g, ''); 
                if (text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
              });
              csv.push(safeCols.join(','));
            }
          });
          downloadCSV(csv.join('\\n'), 'readstats_export.csv');
        }
        
        function exportCoverageToCSV() {
          const table = document.getElementById('dataTable');
          let csv = [];
          const headers = ['Sample', 'Chr', 'Gene/Region', 'Region_Size', 'Pct_Ge_1x', 'Pct_Ge_10x', 'Pct_Ge_20x', 'Pct_Ge_30x'];
          csv.push(headers.join(','));
          const dataRows = table.querySelectorAll('tbody tr');
          dataRows.forEach(row => {
            if (row.style.display !== 'none') {
              const cols = [
                row.getAttribute('data-sample'), row.getAttribute('data-chr'), row.getAttribute('data-gene'),
                row.getAttribute('data-total'), row.getAttribute('data-cov1'), row.getAttribute('data-cov10'),
                row.getAttribute('data-cov20'), row.getAttribute('data-cov30')
              ];
              const safeCols = cols.map(text => {
                text = text.trim();
                if (text.includes(',') || text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
              });
              csv.push(safeCols.join(','));
            }
          });
          downloadCSV(csv.join('\\n'), 'coverage_breadth_export.csv');
        }

        function exportSamtoolsToCSV() {
          const table = document.getElementById('samtoolsTable');
          if (!table) { alert('Coverage table not found.'); return; }
          let csv = [];
          const headers = ['Sample', 'Primary_Mapped', 'Primary_Mapped_Pct', 'Bases_On_Target', 'Mean_Target_Coverage', 'Bases_On_Non_Target', 'Mean_Non_Target_Coverage'];
          csv.push(headers.join(','));
          const dataRows = table.querySelectorAll('tbody tr');
          dataRows.forEach(row => {
            if (row.style.display !== 'none') {
              const cols = [
                row.getAttribute('data-sample'), 
                row.getAttribute('data-pmapped'), row.getAttribute('data-ppct'),
                row.getAttribute('data-tbases'), row.getAttribute('data-tcov'),
                row.getAttribute('data-ntbases'), row.getAttribute('data-ntcov')
              ];
              const safeCols = cols.map(text => {
                if (!text) return "";
                text = text.trim().replace(/,/g, '');
                if (text.includes('"')) { text = '"' + text.replace(/"/g, '""') + '"'; }
                return text;
              });
              csv.push(safeCols.join(','));
            }
          });
          downloadCSV(csv.join('\\n'), 'coverage_stats_export.csv');
        }

        $(document).ready(function() {
          $('#sampleFilter').select2({ placeholder: "Filter samples...", allowClear: true, closeOnSelect: true })
            .on('change', function() { filterTable(); filterReadstatsTable(); filterSamtoolsTable(); }); 
          $('#regionFilter').select2({ placeholder: "Filter regions...", allowClear: true, closeOnSelect: true })
            .on('change', filterTable);
          
          filterTable();
          filterReadstatsTable(); 
          filterSamtoolsTable();
          
          sortDirection[0] = 'desc'; 
          readstatsSortDirection[0] = 'desc';
          samtoolsSortDirection[0] = 'desc';
          sortTable(0);
          sortReadstatsTable(0);
          sortSamtoolsTable(0);
        });

        function sortReadstatsTable(columnIndex) {
          const table = document.getElementById('readstatsTable');
          if (!table) return;
          const tbody = table.tBodies[0];
          const rows = Array.from(tbody.querySelectorAll('tr'));
          readstatsSortDirection[columnIndex] = readstatsSortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
          const isAsc = readstatsSortDirection[columnIndex] === 'asc';
          
          const headers = table.querySelectorAll('th.sortable');
          headers.forEach(h => h.classList.remove('asc', 'desc'));
          if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');
          
          rows.sort((a, b) => {
            const dataAttrMap = { 0: 'sample', 1: 'reads', 2: 'bases', 4: 'minlen', 5: 'maxlen', 6: 'n50', 7: 'gc', 8: 'q20' };
            const dataKey = dataAttrMap[columnIndex];
            if (!dataKey) {
                 const aVal = a.cells[columnIndex].textContent.toLowerCase().trim();
                 const bVal = b.cells[columnIndex].textContent.toLowerCase().trim();
                 return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
            }
            let aVal = a.getAttribute('data-' + dataKey);
            let bVal = b.getAttribute('data-' + dataKey);
            if (columnIndex === 0) return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
            return isAsc ? parseFloat(aVal) - parseFloat(bVal) : parseFloat(bVal) - parseFloat(aVal);
          });
          rows.forEach(row => tbody.appendChild(row));
        }

        function filterReadstatsTable() {
            const selectedSamples = $('#sampleFilter').val() || [];
            const table = document.getElementById('readstatsTable');
            if (!table) return;
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach(row => {
                const sample = row.getAttribute('data-sample');
                const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
                row.style.display = showSample ? '' : 'none';
            });
        }

        function filterSamtoolsTable() {
            const selectedSamples = $('#sampleFilter').val() || [];
            const table = document.getElementById('samtoolsTable');
            if (!table) return;
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach(row => {
                const sample = row.getAttribute('data-sample');
                const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
                row.style.display = showSample ? '' : 'none';
            });
        }

        function sortSamtoolsTable(columnIndex) {
          const table = document.getElementById('samtoolsTable');
          if (!table) return;
          const tbody = table.tBodies[0];
          const rows = Array.from(tbody.querySelectorAll('tr'));
          samtoolsSortDirection[columnIndex] = samtoolsSortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
          const isAsc = samtoolsSortDirection[columnIndex] === 'asc';
          
          const headers = table.querySelectorAll('th.sortable');
          headers.forEach(h => h.classList.remove('asc', 'desc'));
          if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');
          
          rows.sort((a, b) => {
            const dataAttrMap = { 
                0: 'sample', 
                1: 'pmapped', 2: 'ppct',
                3: 'tbases', 4: 'tcov', 
                5: 'ntbases', 6: 'ntcov' 
            };
            const dataKey = dataAttrMap[columnIndex];
            if (columnIndex === 0) {
                 const aVal = a.getAttribute('data-sample');
                 const bVal = b.getAttribute('data-sample');
                 return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
            }
            let aVal = parseFloat(a.getAttribute('data-' + dataKey));
            let bVal = parseFloat(b.getAttribute('data-' + dataKey));
            return isAsc ? aVal - bVal : bVal - aVal;
          });
          rows.forEach(row => tbody.appendChild(row));
        }
        
        function sortTable(columnIndex) {
            const table = document.getElementById('dataTable');
            if (!table) return;
            const tbody = table.tBodies[0];
            const rows = Array.from(tbody.querySelectorAll('tr'));
            sortDirection[columnIndex] = sortDirection[columnIndex] === 'asc' ? 'desc' : 'asc';
            const isAsc = sortDirection[columnIndex] === 'asc';
            
            const headers = table.querySelectorAll('th.sortable');
            headers.forEach(h => h.classList.remove('asc', 'desc'));
            if (headers[columnIndex]) headers[columnIndex].classList.add(isAsc ? 'asc' : 'desc');
            
            rows.sort((a, b) => {
                const dataAttrMap = { 0: 'sample', 1: 'chr', 2: 'gene', 3: 'total', 4: 'cov1', 5: 'cov10', 6: 'cov20', 7: 'cov30' };
                const dataKey = dataAttrMap[columnIndex];
                let aVal = a.getAttribute('data-' + dataKey);
                let bVal = b.getAttribute('data-' + dataKey);
                
                if ([0, 1, 2].includes(columnIndex)) {
                    return isAsc ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
                } else {
                    return isAsc ? parseFloat(aVal) - parseFloat(bVal) : parseFloat(bVal) - parseFloat(aVal);
                }
            });
            rows.forEach(row => tbody.appendChild(row));
        }

        function filterTable() {
            const selectedSamples = $('#sampleFilter').val() || [];
            const selectedRegions = $('#regionFilter').val() || [];
            const table = document.getElementById('dataTable');
            const rows = table.querySelectorAll('tbody tr');
            rows.forEach(row => {
                const sample = row.getAttribute('data-sample');
                const region = row.getAttribute('data-gene');
                const showSample = selectedSamples.length === 0 || selectedSamples.includes(sample);
                const showRegion = selectedRegions.length === 0 || selectedRegions.includes(region);
                row.style.display = (showSample && showRegion) ? '' : 'none';
            });
        }
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
              <td style="padding: 3px 10px 3px 0; font-weight: 500; border: none; background: inherit; font-family: 'Courier New', monospace; font-weight: 600; color: white; white-space: nowrap; font-size: 0.8em; width: 200px;">{display_key}:</td>
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

def render_stats_cards(readstats_data, samples_data, genes, region_totals, ref_stats):
    """Render the top statistics cards"""
    total_bases = sum(stats.get('bases', 0) for stats in readstats_data.values()) if readstats_data else 0
    total_reads = sum(stats.get('reads', 0) for stats in readstats_data.values()) if readstats_data else 0
    total_bed_size = sum(region_totals.values())
    
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
          <h3>Samples</h3>
          <div class="value">{sample_count}</div>
        </div>
        <div class="stat-card">
          <h3>Total Reads</h3>
          <div class="value">{reads_fmt}</div>
        </div>
        <div class="stat-card">
          <h3>Total Bases</h3>
          <div class="value">{bases_fmt}</div>
        </div>
        <div class="stat-card">
          <h3>Ref Size</h3>
          <div class="value">{ref_bases_fmt}</div>
        </div>
        <div class="stat-card">
          <h3>Ref Contigs</h3>
          <div class="value">{ref_contigs}</div>
        </div>
        <div class="stat-card">
          <h3>BED Genes/Regions</h3>
          <div class="value">{genes_count}</div>
        </div>
        <div class="stat-card">
          <h3>Total BED size</h3>
          <div class="value">{bed_fmt}</div>
        </div>
      </div>
    """
    return html

def render_readstats_table(readstats_data):
    """Render the Read Statistics table"""
    if not readstats_data:
        return ""
    
    html = """
      <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
        <h3 style="margin: 0; font-size: 1.1em; color: #374151;">Read Statistics</h3>
        <button class="export-btn" onclick="exportReadstatsToCSV()">Export CSV</button>
      </div>
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
        
    html += """
          </tbody>
        </table>
      </div>
      
    """
    return html

def render_samtools_table(samtools_data):
    """Render the Samtools Coverage Statistics table"""
    if not samtools_data:
        return ""

    html = """
      <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
        <h3 style="margin: 0; font-size: 1.1em; color: #374151;">Coverage</h3>
        <button class="export-btn" onclick="exportSamtoolsToCSV()">Export CSV</button>
      </div>
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
    """
    return html

def render_coverage_table(samples_data, genes):
    """Render the Coverage Statistics table"""
    html = """
      <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 15px;">
        <h3 style="margin: 0; font-size: 1.1em; color: #374151;">Breadth of Coverage</h3>
        <button class="export-btn" onclick="exportCoverageToCSV()">Export CSV</button>
      </div>
      <div class="table-container">
        <table id="dataTable">
          <thead>
            <tr>
              <th rowspan="2" class="sample-col sortable" onclick="sortTable(0)">Sample</th>
              <th rowspan="2" class="sortable" onclick="sortTable(1)">Chr</th>
              <th rowspan="2" class="sortable" onclick="sortTable(2)">Gene/Region</th>
              <th rowspan="2" class="sortable" onclick="sortTable(3)" style="text-align: right;">Region size</th>
              <th colspan="4" style="text-align: center; border-bottom: 1px solid #cbd5e1;">Percentage of region with at least X coverage</th>
            </tr>
            <tr>
              <th class="sortable" onclick="sortTable(4)" style="text-align: center;">≥1x (%)</th>
              <th class="sortable" onclick="sortTable(5)" style="text-align: center;">≥10x (%)</th>
              <th class="sortable" onclick="sortTable(6)" style="text-align: center;">≥20x (%)</th>
              <th class="sortable" onclick="sortTable(7)" style="text-align: center;">≥30x (%)</th>
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
              <td>{location['chr']}</td>
              <td><strong>{gene}</strong></td>
              <td style="text-align: right;">{cov_stats['total']:,}</td>
              <td class="coverage-cell">
                <span class="pct-bar" style="width: {cov_stats['pct_1x']*1.5}px;">{cov_stats['pct_1x']:.1f}%</span>
              </td>
              <td class="coverage-cell">
                <span class="pct-bar" style="width: {cov_stats['pct_10x']*1.5}px;">{cov_stats['pct_10x']:.1f}%</span>
              </td>
              <td class="coverage-cell">
                <span class="pct-bar" style="width: {cov_stats['pct_20x']*1.5}px;">{cov_stats['pct_20x']:.1f}%</span>
              </td>
              <td class="coverage-cell">
                <span class="pct-bar" style="width: {cov_stats['pct_30x']*1.5}px;">{cov_stats['pct_30x']:.1f}%</span>
              </td>
            </tr>
            """
    
    html += """
          </tbody>
        </table>
      </div>
    """
    return html

def generate_html_report(samples_data, readstats_data, run_info, wf_info, ref_stats, samtools_stats, output_file):
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
    
    sample_names = sorted(list(all_sample_names))
    sample_options = "".join([f'<option value="{name.lower()}">{name}</option>' for name in sample_names])

    # Render Components
    css_block = get_css()
    js_block = get_js()
    
    run_info_block = render_details_block("Sequencing run details", run_info, add_top_border=True)
    wf_info_block = render_details_block("Workflow details", wf_info, add_top_border=False)
    
    stats_cards = render_stats_cards(readstats_data, samples_data, genes, region_totals, ref_stats)
    stats_cards = render_stats_cards(readstats_data, samples_data, genes, region_totals, ref_stats)
    readstats_table = render_readstats_table(readstats_data)
    samtools_table = render_samtools_table(samtools_stats)
    coverage_table = render_coverage_table(samples_data, genes)
    
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
      <h2>NXF-ALIGNMENT Report</h2>
      <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
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
      {coverage_table}
    </div>
  </div>
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
    parser.add_argument('-o', '--output', required=True, help='Output HTML file')
    
    args = parser.parse_args()
    
    samples_data = {}
    readstats_data = {}
    samtools_stats = {}
    run_info = [] 
    wf_info = []
    ref_stats = []

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
        print(f"Processing stats {readstats_file} (Sample: {sample_name})...")
        readstats_data[sample_name] = parse_readstats_file(readstats_file)

    # Process bedcov files
    # First, collect all potential samples
    bedcov_map = {Path(f).name.replace('.batched.reads.bedcov.tsv', '').replace('.reads.bedcov.tsv', '').replace('.bedcov.tsv', ''): f for f in args.bedcov}
    bedcov_comp_map = {Path(f).name.replace('.batched.reads.bedcov.compl.tsv', '').replace('.reads.bedcov.compl.tsv', '').replace('.bedcov.compl.tsv', ''): f for f in args.bedcov_compl}
    
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
        if name not in samtools_stats: samtools_stats[name] = {}
        samtools_stats[name]['target'] = parse_bedcov_file(f)
        
    for f in args.bedcov_compl:
        name = Path(f).name
        for suffix in ['.batched.reads.bedcov.compl.tsv', '.reads.bedcov.compl.tsv', '.bedcov.compl.tsv']:
             if name.endswith(suffix):
                 name = name[:-len(suffix)]
                 break
        if name not in samtools_stats: samtools_stats[name] = {}
        samtools_stats[name]['non-target'] = parse_bedcov_file(f)

    for f in args.flagstat:
        name = Path(f).name
        # Assuming filename format like sample.flagstat.json
        for suffix in ['.flagstat.json']:
             if name.endswith(suffix):
                 name = name[:-len(suffix)]
                 break
        if name not in samtools_stats: samtools_stats[name] = {}
        samtools_stats[name]['flagstat'] = parse_flagstat_file(f)

    generate_html_report(samples_data, readstats_data, run_info, wf_info, ref_stats, samtools_stats, args.output)

if __name__ == "__main__":
    main()