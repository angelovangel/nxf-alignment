#!/usr/bin/env python3
"""
Generate an HTML report from snpEff stats CSV files.
Usage: python make-variants-report.py csv1.stats.csv [csv2.stats.csv ...] -o output.html
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

def parse_snpeff_csv(filepath):
    """Parse a snpEff stats.csv file"""
    data = {
        'summary': {},
        'change_rate_by_chr': [],
        'variant_types': [],
        'effects_impact': [],
        'functional_class': [],
        'count_by_effects': [],
        'count_by_region': [],
        'quality_dist': { 'values': [], 'counts': [] },
        'indel_lengths': { 'values': [], 'counts': [] },
        'base_changes': [], # Matrix
        'tstv': {},
        'zygosity': {}
    }
    
    current_section = None
    
    with open(filepath, 'r') as f:
        reader = csv.reader(f, skipinitialspace=True)
        for row in reader:
            if not row:
                continue
            
            line = ",".join(row).strip()
            if line.startswith("# Summary table"):
                current_section = "summary"
                continue
            elif line.startswith("# Change rate by chromosome"):
                current_section = "change_rate"
                continue
            elif line.startswith("# Variantss by type"):
                current_section = "variant_types"
                continue
            elif line.startswith("# Effects by impact"):
                current_section = "effects_impact"
                continue
            elif line.startswith("# Effects by functional class"):
                current_section = "functional_class"
                continue
            elif line.startswith("# Count by effects"):
                current_section = "count_by_effects"
                continue
            elif line.startswith("# Count by genomic region"):
                current_section = "count_by_region"
                continue
            elif line.startswith("# Quality"):
                current_section = "quality"
                continue
            elif line.startswith("# InDel lengths"):
                current_section = "indel_lengths"
                continue
            elif line.startswith("# Base changes"):
                current_section = "base_changes"
                continue
            elif line.startswith("# Ts/Tv summary"):
                current_section = "tstv"
                continue
            elif line.startswith("# Hom/Het table"):
                current_section = "zygosity"
                continue
            elif line.startswith("# Codon change table") or line.startswith("# Amino acid change table"):
                current_section = "ignore"
                continue

            if current_section == "ignore":
                continue

            # Section parsing
            if current_section == "summary":
                if len(row) >= 2 and row[0].strip() not in ["Name", ""]:
                    data['summary'][row[0].strip()] = row[1].strip()
                elif len(row) == 1 and "Missense_Silent_ratio" in row[0]:
                    parts = row[0].split(",")
                    if len(parts) >= 2:
                        try: data['summary']['Missense_Silent_ratio'] = parts[1].strip()
                        except: pass

            elif current_section == "change_rate":
                if len(row) >= 4 and row[0].strip() not in ["Chromosome", ""]:
                    try:
                        data['change_rate_by_chr'].append({
                            'chr': row[0].strip(),
                            'length': int(row[1].strip()),
                            'changes': int(row[2].strip()),
                            'rate': float(row[3].strip())
                        })
                    except ValueError: pass

            elif current_section == "variant_types":
                if len(row) >= 2 and row[0].strip() not in ["Type", ""]:
                    try:
                        data['variant_types'].append({
                            'type': row[0].strip(),
                            'count': int(row[1].strip()),
                            'percent': row[2].strip() if len(row) > 2 else "0%"
                        })
                    except ValueError: pass

            elif current_section == "effects_impact":
                if len(row) >= 2 and row[0].strip() not in ["Type", ""]:
                    try:
                        data['effects_impact'].append({
                            'type': row[0].strip(),
                            'count': int(row[1].strip()),
                            'percent': row[2].strip() if len(row) > 2 else "0%"
                        })
                    except ValueError: pass

            elif current_section == "functional_class":
                if "Missense_Silent_ratio" in row[0]:
                    try:
                        data['summary']['Missense_Silent_ratio'] = row[1].strip()
                    except (IndexError, ValueError): pass
                elif len(row) >= 2 and row[0].strip() not in ["Type", ""]:
                    try:
                        data['functional_class'].append({
                            'type': row[0].strip(),
                            'count': int(row[1].strip()),
                            'percent': row[2].strip() if len(row) > 2 else "0%"
                        })
                    except ValueError: pass

            elif current_section == "count_by_effects":
                if len(row) >= 2 and row[0].strip() not in ["Type", ""]:
                    try:
                        data['count_by_effects'].append({
                            'type': row[0].strip(),
                            'count': int(row[1].strip()),
                            'percent': row[2].strip() if len(row) > 2 else "0%"
                        })
                    except ValueError: pass

            elif current_section == "count_by_region":
                if len(row) >= 2 and row[0].strip() not in ["Type", ""]:
                    try:
                        data['count_by_region'].append({
                            'type': row[0].strip(),
                            'count': int(row[1].strip()),
                            'percent': row[2].strip() if len(row) > 2 else "0%"
                        })
                    except ValueError: pass

            elif current_section == "quality":
                if row[0].strip() == "Values":
                    data['quality_dist']['values'] = [v.strip() for v in row[1:]]
                elif row[0].strip() == "Count":
                    data['quality_dist']['counts'] = [int(v.strip()) for v in row[1:]]

            elif current_section == "indel_lengths":
                if row[0].strip() == "Values":
                    data['indel_lengths']['values'] = [v.strip() for v in row[1:]]
                elif row[0].strip() == "Count":
                    data['indel_lengths']['counts'] = [int(v.strip()) for v in row[1:]]

            elif current_section == "base_changes":
                if row[0].strip() == "base":
                    data['base_header'] = [r.strip() for r in row[1:] if r.strip()]
                elif len(row) > 1 and row[0].strip() in ["A", "C", "G", "T"]:
                    try:
                        data['base_changes'].append({
                            'ref': row[0].strip(),
                            'counts': [int(v.strip()) for v in row[1:] if v.strip()]
                        })
                    except ValueError: pass

            elif current_section == "tstv":
                if len(row) >= 2:
                    data['tstv'][row[0].strip()] = row[1].strip()

            elif current_section == "zygosity":
                if len(row) >= 2 and row[0].strip() not in ["Sample_names", ""]:
                    data['zygosity'][row[0].strip()] = row[1].strip()

    return data

def get_css():
    css_path = Path(__file__).parent / 'report' / 'assets' / 'report.css'
    if css_path.exists():
        with open(css_path, 'r') as f:
            css_content = f.read()
    else:
        css_content = "body { font-family: sans-serif; }" # Minimal fallback

    return f"""
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;600;700&display=swap" rel="stylesheet">
    <style>
{css_content}
    .grid-2 {{ display: grid; grid-template-columns: 1fr 1fr; gap: 20px; }}
    .impact-HIGH {{ color: #dc2626; font-weight: bold; font-family: 'Inter', sans-serif; }}
    .impact-MODERATE {{ color: #ea580c; font-family: 'Inter', sans-serif; }}
    .impact-LOW {{ color: #16a34a; font-family: 'Inter', sans-serif; }}
    .impact-MODIFIER {{ color: #64748b; font-family: 'Inter', sans-serif; }}
    .base-matrix {{ border-collapse: collapse; margin: 10px 0; font-family: monospace; }}
    .base-matrix td, .base-matrix th {{ padding: 8px; border: 1px solid #e2e8f0; text-align: right; }}
    .base-matrix th {{ background: #f8fafc; }}
    .base-matrix .ref-col {{ font-weight: bold; background: #f8fafc; text-align: center; }}
    .dist-bar {{ display: inline-block; height: 10px; background: #374151; border-radius: 2px; }}
    </style>
    """

def get_js():
    js_path = Path(__file__).parent / 'report' / 'assets' / 'report.js'
    if js_path.exists():
        with open(js_path, 'r') as f:
            js_content = f.read()
    else:
        js_content = ""
    return f"""
    <script src="https://ajax.googleapis.com/ajax/libs/jquery/3.5.1/jquery.min.js"></script>
    <script>{js_content}</script>
    """

def render_summary_cards(all_data):
    total_variants = sum(int(d['summary'].get('Number_of_variants_processed', '0')) for d in all_data.values())
    total_effects = sum(int(d['summary'].get('Number_of_effects', '0')) for d in all_data.values())
    
    html = f"""
    <div class="stats">
        <div class="stat-card"><h3>Samples</h3><div class="value">{len(all_data)}</div></div>
        <div class="stat-card"><h3>Total Variants</h3><div class="value">{format_si(total_variants)}</div></div>
        <div class="stat-card"><h3>Total Effects</h3><div class="value">{format_si(total_effects)}</div></div>
    </div>
    """
    return html

def render_main_table(all_data):
    html = """
    <details class="collapsible-section" open>
        <summary><h3>Sample Summary</h3></summary>
        <div class="table-container">
            <table>
                <thead>
                    <tr>
                        <th>Sample</th>
                        <th>Variants</th>
                        <th>SNP</th>
                        <th>INS</th>
                        <th>DEL</th>
                        <th>Ts/Tv</th>
                        <th>Het</th>
                        <th>Hom</th>
                        <th>Missense</th>
                        <th>Nonsense</th>
                        <th>Silent</th>
                        <th>Miss/Sil</th>
                    </tr>
                </thead>
                <tbody>
    """
    for sample, data in sorted(all_data.items()):
        v_types = {t['type']: t['count'] for t in data['variant_types']}
        f_class = {t['type']: t['count'] for t in data['functional_class']}
        html += f"""
        <tr>
            <td><strong>{sample}</strong></td>
            <td>{data['summary'].get('Number_of_variants_processed', '0')}</td>
            <td>{v_types.get('SNP', 0):,}</td>
            <td>{v_types.get('INS', 0):,}</td>
            <td>{v_types.get('DEL', 0):,}</td>
            <td>{data['tstv'].get('Ts_Tv_ratio', 'N/A')}</td>
            <td>{data['zygosity'].get('Het', '0')}</td>
            <td>{data['zygosity'].get('Hom', '0')}</td>
            <td>{f_class.get('MISSENSE', 0):,}</td>
            <td>{f_class.get('NONSENSE', 0):,}</td>
            <td>{f_class.get('SILENT', 0):,}</td>
            <td>{data['summary'].get('Missense_Silent_ratio', 'N/A')}</td>
        </tr>"""
    html += "</tbody></table></div></details>"
    return html

def render_detailed_section(sample, data):
    # Effects Table
    effects_rows = "".join([f"<tr><td>{e['type']}</td><td>{e['count']:,}</td><td>{e['percent']}</td></tr>" for e in data['count_by_effects']])
    regions_rows = "".join([f"<tr><td>{e['type']}</td><td>{e['count']:,}</td><td>{e['percent']}</td></tr>" for e in data['count_by_region']])
    impacts_rows = "".join([f"<tr><td class='impact-{e['type']}'>{e['type']}</td><td>{e['count']:,}</td><td>{e['percent']}</td></tr>" for e in data['effects_impact']])
    
    # Base Changes Matrix
    base_html = "<table class='base-matrix'><thead><tr><th>Ref \\ Alt</th>"
    base_html += "".join([f"<th>{b}</th>" for b in data.get('base_header', [])])
    base_html += "</tr></thead><tbody>"
    for row in data['base_changes']:
        base_html += f"<tr><td class='ref-col'>{row['ref']}</td>"
        base_html += "".join([f"<td>{c:,}</td>" for c in row['counts']])
        base_html += "</tr>"
    base_html += "</tbody></table>"

    safe_sample = sample.replace("-", "_").replace(".", "_").replace(" ", "_")

    html = f"""
    <details class="collapsible-section" open>
        <summary><h3>Details for Sample: {sample}</h3></summary>
        <div style="padding: 20px;">
            <div class="grid-2">
                <div class="collapsible-section">
                    <summary style="background: #374151; color: white;">Quality & InDel Distributions</summary>
                    <div style="padding: 15px; background: #f8fafc; border-radius: 8px;">
                        <div style="height: 200px; margin-bottom: 20px;">
                            <canvas id="qualityChart_{safe_sample}"></canvas>
                        </div>
                        <div style="height: 200px;">
                            <canvas id="indelChart_{safe_sample}"></canvas>
                        </div>
                    </div>
                </div>
                <div class="collapsible-section">
                    <summary style="background: #374151; color: white;">Base Changes (Ref vs Alt)</summary>
                    <div style="padding: 10px;">{base_html}</div>
                </div>
            </div>

            <div class="grid-2" style="margin-top: 20px;">
                <div class="collapsible-section">
                    <summary style="background: #374151; color: white;">Effects by Type</summary>
                    <div class="table-container"><table><thead><tr><th>Effect</th><th>Count</th><th>%</th></tr></thead><tbody>{effects_rows}</tbody></table></div>
                </div>
                <div class="collapsible-section">
                    <summary style="background: #374151; color: white;">Effects by Region</summary>
                    <div class="table-container"><table><thead><tr><th>Region</th><th>Count</th><th>%</th></tr></thead><tbody>{regions_rows}</tbody></table></div>
                </div>
            </div>
            
            <div class="grid-2" style="margin-top: 20px;">
                <div class="collapsible-section">
                    <summary style="background: #374151; color: white;">Effects by Impact</summary>
                    <div class="table-container"><table><thead><tr><th>Impact</th><th>Count</th><th>%</th></tr></thead><tbody>{impacts_rows}</tbody></table></div>
                </div>
                <div class="collapsible-section" style="margin-bottom: 0;">
                    <summary style="background: #374151; color: white;">Change Rate by Chromosome</summary>
                    <div class="table-container">
                        <table><thead><tr><th>Chr</th><th>Length</th><th>Changes</th><th>Rate (1 / bp)</th></tr></thead>
                        <tbody>
                            {"".join([f"<tr><td>{c['chr']}</td><td>{c['length']:,}</td><td>{c['changes']:,}</td><td>{c['rate']:,}</td></tr>" for c in data['change_rate_by_chr']])}
                        </tbody></table>
                    </div>
                </div>
            </div>
        </div>
    </details>

    <script>
    (function() {{
        initVariantCharts(
            '{safe_sample}', 
            {json.dumps(data['quality_dist']['values'])}, 
            {json.dumps(data['quality_dist']['counts'])},
            {json.dumps(data['indel_lengths']['values'])}, 
            {json.dumps(data['indel_lengths']['counts'])}
        );
    }})();
    </script>
    """
    return html

def main():
    parser = argparse.ArgumentParser(description="Generate HTML report from snpEff stats CSV")
    parser.add_argument("csv_files", nargs="+", help="One or more snpEff stats CSV files")
    parser.add_argument("-o", "--output", default="variants_report.html", help="Output HTML file")
    args = parser.parse_args()
    
    all_data = {}
    for csv_file in args.csv_files:
        path = Path(csv_file)
        sample = path.name.replace(".stats.csv", "").replace(".csv", "")
        all_data[sample] = parse_snpeff_csv(csv_file)
    
    if not all_data:
        print("No data found.")
        return

    css_block = get_css()
    js_block = get_js()
    stats_cards = render_summary_cards(all_data)
    summary_table = render_main_table(all_data)
    
    detailed_sections = "".join([render_detailed_section(s, d) for s, d in sorted(all_data.items())])
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Variant Annotation Report</title>
    {css_block}
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    {js_block}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-main"><h2>Variant Annotation Report</h2></div>
            <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        <div class="content">
            {stats_cards}
            {summary_table}
            {detailed_sections}
        </div>
    </div>
</body>
</html>
"""
    with open(args.output, "w") as f: f.write(html_content)
    print(f"Report generated: {args.output}")

if __name__ == "__main__":
    main()
