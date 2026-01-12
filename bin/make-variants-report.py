#!/usr/bin/env python3
"""
Generate an HTML report from snpEff stats CSV files.
Usage: python make-variants-report.py csv1.stats.csv [csv2.stats.csv ...] -o output.html
"""

import argparse
from pathlib import Path
from datetime import datetime
import csv

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

def fmt_n(val):
    """Format large numbers with commas"""
    try:
        return f"{int(str(val).replace(',', '')):,.0f}"
    except (ValueError, TypeError):
        return str(val)



def parse_snpeff_csv(filepath):
    """Parse a snpEff stats.csv file"""
    data = {
        'summary': {},
        'variant_types': [],
        'effects_impact': [],
        'functional_class': [],
        'count_by_effects': [],
        'count_by_region': [],
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
            elif line.startswith("# Ts/Tv summary"):
                current_section = "tstv"
                continue
            elif line.startswith("# Hom/Het table"):
                current_section = "zygosity"
                continue
            elif any(line.startswith(p) for p in [
                "# Change rate by chromosome", 
                "# Quality", 
                "# InDel lengths", 
                "# Base changes", 
                "# Codon change table", 
                "# Amino acid change table"
            ]):
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

def render_header_details(all_data):
    """Render details block for the header using summary data from the first sample"""
    if not all_data:
        return ""
    
    # Get summary from the first sample
    first_sample = list(all_data.values())[0]
    summary = first_sample.get('summary', {})
    
    # Fields to display
    display_fields = ["SnpEff_version", "Genome", "Date", "Command_line_arguments"]
    
    html = """
    <details style="padding-top: 5px; margin-top: 5px; border-top: 1px solid rgba(255, 255, 255, 0.2); text-align: left;">
      <summary style="font-size: 0.8em; cursor: pointer; color: white;">SnpEff Details</summary>
      <div style="padding-top: 5px; text-align: left; width: 100%;">
        <table style="width: 100%; border-collapse: collapse; margin-top: 5px; color: white; border: none; background: inherit;">
          <tbody>
    """
    
    for key in display_fields:
        if key in summary:
            display_key = key.replace('_', ' ').title()
            value = summary[key]
            html += f"""
            <tr>
              <td style="padding: 3px 10px 3px 0; font-weight: 500; border: none; background: inherit; font-family: 'Courier New', monospace; color: white; white-space: nowrap; font-size: 0.8em; width: 180px;">{display_key}:</td>
              <td style="padding: 3px 0; border: none; background: inherit; color: rgba(255, 255, 255, 0.8); font-size: 0.8em;">{value}</td>
            </tr>
            """
            
    html += """
          </tbody>
        </table>
      </div>
    </details>
    """
    return html

def render_aggregate_tables(all_data):
    """Render consolidated tables for Effects, Regions, and Impacts"""
    sections = {
        'count_by_effects': 'Effects by Type',
        'count_by_region': 'Effects by Region',
        'effects_impact': 'Effects by Impact'
    }
    
    
    html = ""
    
    # Abbreviations for Effects by Type
    effect_abbr = {
        '3_prime_UTR_variant': "3' UTR",
        '5_prime_UTR_premature_start_codon_gain_variant': "5' UTR Stop+",
        '5_prime_UTR_variant': "5' UTR",
        'disruptive_inframe_deletion': 'Dis. Inframe Del',
        'disruptive_inframe_insertion': 'Dis. Inframe Ins',
        'downstream_gene_variant': 'Downstream',
        'intergenic_region': 'Intergenic',
        'intron_variant': 'Intron',
        'missense_variant': 'Missense',
        'non_coding_transcript_exon_variant': 'NC Exon',
        'sequence_feature': 'Seq Feat',
        'splice_acceptor_variant': 'Splice Acc',
        'splice_donor_variant': 'Splice Don',
        'splice_region_variant': 'Splice Reg',
        'stop_gained': 'Stop Gained',
        'stop_lost': 'Stop Lost',
        'stop_retained_variant': 'Stop Retained',
        'structural_interaction_variant': 'Struct Int',
        'synonymous_variant': 'Synonymous',
        'upstream_gene_variant': 'Upstream',
        'start_lost': 'Start Lost',
        'initiator_codon_variant': 'Init Codon'
    }

    for key, title in sections.items():
        # Collect all unique columns (types)
        all_types = set()
        for sample, data in all_data.items():
            for item in data.get(key, []):
                all_types.add(item['type'])
        sorted_types = sorted(list(all_types))
        
        if not sorted_types:
            continue
            
        header_row = ""
        for t in sorted_types:
            if key == 'count_by_effects':
                display_text = effect_abbr.get(t, t.replace('_', ' ').title())
                header_row += f"<th title='{t}'>{display_text}</th>"
            else:
                display_text = t.replace('_', ' ').title()
                header_row += f"<th title='{t}'>{display_text}</th>"

        html += f"""
        <details class="collapsible-section" open>
            <summary><h3>{title}</h3></summary>
            <div class="table-container">
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            {header_row}
                        </tr>
                    </thead>
                    <tbody>
        """
        
        for sample, data in sorted(all_data.items()):
            # Create a lookup for this sample's counts
            counts = {item['type']: item['count'] for item in data.get(key, [])}
            
            html += f"<tr><td><strong>{sample}</strong></td>"
            for t in sorted_types:
                val = counts.get(t, 0)
                # Calculate max for this column across all samples for heatmap?
                # User asked for heatmap on percentages only previously. 
                # For this consolidated table, showing raw counts is standard.
                html += f"<td>{fmt_n(val)}</td>"
            html += "</tr>"
            
        html += """
                    </tbody>
                </table>
            </div>
        </details>
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
            <td>{fmt_n(data['summary'].get('Number_of_variants_processed', 0))}</td>
            <td>{fmt_n(v_types.get('SNP', 0))}</td>
            <td>{fmt_n(v_types.get('INS', 0))}</td>
            <td>{fmt_n(v_types.get('DEL', 0))}</td>
            <td>{data['tstv'].get('Ts_Tv_ratio', 'N/A')}</td>
            <td>{fmt_n(data['zygosity'].get('Het', 0))}</td>
            <td>{fmt_n(data['zygosity'].get('Hom', 0))}</td>
            <td>{fmt_n(f_class.get('MISSENSE', 0))}</td>
            <td>{fmt_n(f_class.get('NONSENSE', 0))}</td>
            <td>{fmt_n(f_class.get('SILENT', 0))}</td>
            <td>{data['summary'].get('Missense_Silent_ratio', 'N/A')}</td>
        </tr>"""
    html += "</tbody></table></div></details>"
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
    header_details = render_header_details(all_data)
    summary_table = render_main_table(all_data)
    aggregate_tables = render_aggregate_tables(all_data)
    
    
    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Variant Annotation Report</title>
    {css_block}

    {js_block}
</head>
<body>
    <div class="container">
        <div class="header">
            <div class="header-main"><h2>Variant Annotation Report</h2></div>
            <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            {header_details}
        </div>
        <div class="content">
            {stats_cards}
            {summary_table}
            {aggregate_tables}
        </div>
    </div>
</body>
</html>
"""
    with open(args.output, "w") as f: f.write(html_content)
    print(f"Report generated: {args.output}")

if __name__ == "__main__":
    main()
