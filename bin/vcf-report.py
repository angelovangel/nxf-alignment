#!/usr/bin/env python3
"""
vcf_report.py — VCF + PharmCAT JSON → HTML Clinical Report
Usage:
    python vcf_report.py --vcf sample.vcf.gz --report pharmcat_report.json [options]

Optional inputs:
    --out     report.html  (default: <sample>_report.html)
"""

import argparse
import gzip
import json
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path


# ─── VCF PARSING ────────────────────────────────────────────────────────────

def open_vcf(path: str):
    p = Path(path)
    if p.suffix in ('.gz', '.bgz') or path.endswith('.vcf.gz'):
        return gzip.open(path, 'rt', encoding='utf-8')
    return open(path, 'r', encoding='utf-8')


def parse_ann(ann_str: str) -> dict:
    """Parse SnpEff ANN= field, return dict from first transcript annotation."""
    # ANN=ALT|consequence|impact|gene|gene_id|feature_type|feature_id|biotype|rank|hgvsc|hgvsp|...
    fields = ann_str.split('|')
    labels = [
        'allele', 'consequence', 'impact', 'gene', 'gene_id',
        'feature_type', 'feature_id', 'biotype', 'rank',
        'hgvsc', 'hgvsp', 'cdna_pos', 'cds_pos', 'aa_pos', 'distance', 'warnings'
    ]
    return {labels[i]: fields[i] if i < len(fields) else '' for i in range(len(labels))}


def parse_info(info_str: str) -> dict:
    result = {}
    for part in info_str.split(';'):
        if '=' in part:
            k, v = part.split('=', 1)
            result[k] = v
        else:
            result[part] = True
    return result


def gt_zygosity(gt: str, ref: str = None) -> str:
    if not gt or gt == './.': return 'REF'
    parts = gt.replace('|', '/').split('/')
    
    def is_ref_allele(p):
        if p == '.': return True
        if p.isdigit(): return p == '0'
        return p == ref if ref else False

    if all(is_ref_allele(p) for p in parts):
        return 'REF'
    
    called = [p for p in parts if p != '.']
    if not called: return 'REF'
    if len(set(called)) == 1 and len(called) > 1:
        return 'HOM'
    return 'HET'


def compute_sample_af(fmt: dict) -> float:
    """Compute alt allele fraction from AD field if FORMAT AF absent."""
    if fmt.get('AF') and fmt['AF'] != '.':
        try:
            return float(fmt['AF'].split(',')[0])
        except ValueError:
            pass
    if fmt.get('AD') and fmt['AD'] != '.':
        parts = fmt['AD'].split(',')
        try:
            ref_ad = int(parts[0])
            alt_ad = sum(int(x) for x in parts[1:] if x.isdigit())
            total = ref_ad + alt_ad
            return alt_ad / total if total > 0 else 0.0
        except (ValueError, IndexError):
            pass
    return 0.0


IMPACT_ORDER = {'HIGH': 0, 'MODERATE': 1, 'LOW': 2, 'MODIFIER': 3}

def acmg_classify(consequence: str, impact: str, info_af: float) -> str:
    """Simplified ACMG-style classification."""
    high_impact = {'stop_gained', 'frameshift_variant', 'splice_acceptor_variant',
                   'splice_donor_variant', 'start_lost', 'stop_lost'}
    moderate_impact = {'missense_variant', 'protein_altering_variant',
                       'splice_region_variant', 'inframe_insertion', 'inframe_deletion'}

    csqs = set(consequence.split('&'))
    if csqs & high_impact:
        return 'Pathogenic' if info_af < 0.001 else 'Likely Pathogenic'
    if csqs & moderate_impact:
        if info_af < 0.001:
            return 'Likely Pathogenic'
        if info_af < 0.01:
            return 'VUS'
        if info_af < 0.05:
            return 'Likely Benign'
        return 'Benign'
    if 'synonymous_variant' in csqs:
        return 'Benign'
    if info_af > 0.05:
        return 'Benign'
    if info_af > 0.01:
        return 'Likely Benign'
    return 'VUS'


def parse_vcf(vcf_path: str) -> tuple[str, dict[tuple[str, int], dict], dict]:
    """Returns (sample_name, variant_map, qc_stats)."""
    variant_map = {}
    sample_name = 'SAMPLE'
    meta_lines = []
    genome_build = 'Unknown'

    with open_vcf(vcf_path) as fh:
        for raw in fh:
            line = raw.rstrip('\n')
            if line.startswith('##'):
                meta_lines.append(line)
                if 'reference' in line.lower() or 'genome' in line.lower():
                    if 'GRCh38' in line or 'hg38' in line:
                        genome_build = 'GRCh38/hg38'
                    elif 'GRCh37' in line or 'hg19' in line:
                        genome_build = 'GRCh37/hg19'
                continue

            if line.startswith('#CHROM'):
                cols = line.split('\t')
                if len(cols) > 9:
                    sample_name = cols[9]
                continue

            if not line.strip():
                continue

            cols = line.split('\t')
            if len(cols) < 8:
                continue

            chrom, pos, vid, ref, alt, qual, flt, info_str = cols[:8]
            fmt_str = cols[8] if len(cols) > 8 else 'GT'
            smp_str = cols[9] if len(cols) > 9 else '.'

            info = parse_info(info_str)
            fmt_keys = fmt_str.split(':')
            smp_vals = smp_str.split(':')
            fmt = {fmt_keys[i]: smp_vals[i] if i < len(smp_vals) else '.' for i in range(len(fmt_keys))}

            gt = fmt.get('GT', './.')
            zyg = gt_zygosity(gt, ref)
            sample_af = compute_sample_af(fmt)

            # Rescue 0/0 calls with alt read support (Clair3/bcftools style)
            if zyg == 'REF' and sample_af > 0.05:
                zyg = 'HET'
            if zyg == 'REF':
                continue

            # Parse ANN field (SnpEff)
            ann = {}
            raw_ann = info.get('ANN') or info.get('CSQ', '')
            if raw_ann:
                # Multiple transcripts separated by comma; pick highest impact
                transcripts = raw_ann.split(',')
                best = None
                best_order = 99
                for t in transcripts:
                    parsed = parse_ann(t)
                    order = IMPACT_ORDER.get(parsed.get('impact', ''), 99)
                    if order < best_order:
                        best_order = order
                        best = parsed
                ann = best or {}

            gene = (ann.get('gene') or info.get('GENE') or info.get('PX') or '').strip()
            consequence = ann.get('consequence', 'unknown')
            impact = ann.get('impact', '')
            hgvsc = ann.get('hgvsc', '')
            hgvsp = ann.get('hgvsp', '')
            transcript = ann.get('feature_id', '')

            # Allele frequency: prefer gnomAD from INFO, else sample AF
            gnomad_af = 0.0
            for af_key in ('AF', 'gnomAD_AF', 'AF_popmax', 'MAX_AF'):
                if info.get(af_key):
                    try:
                        gnomad_af = float(str(info[af_key]).split(',')[0])
                        break
                    except ValueError:
                        pass
            if gnomad_af == 0.0:
                gnomad_af = sample_af

            var_type = 'SNV' if len(ref) == 1 and len(alt.split(',')[0]) == 1 else 'INDEL'

            classification = acmg_classify(consequence, impact, gnomad_af)

            dp = fmt.get('DP', '.')
            gq = fmt.get('GQ', '.')
            ps = fmt.get('PS', '')

            variant_entry = {
                'chr': chrom.lstrip('chr') if chrom.startswith('chr') else chrom,
                'chrom': chrom,
                'pos': int(pos),
                'id': vid if vid != '.' else '',
                'ref': ref,
                'alt': alt,
                'qual': qual,
                'filter': flt,
                'gene': gene,
                'consequence': consequence,
                'impact': impact,
                'hgvsc': hgvsc,
                'hgvsp': hgvsp,
                'transcript': transcript,
                'zygosity': zyg,
                'var_type': var_type,
                'gnomad_af': gnomad_af,
                'sample_af': sample_af,
                'classification': classification,
                'dp': dp,
                'gq': gq,
                'gt': gt,
                'ps': ps,
            }
            variant_map[(chrom, pos)] = variant_entry

    # QC stats
    variants = list(variant_map.values())
    total = len(variants)
    qc = {
        'total': total,
        'sample': sample_name,
        'genome_build': genome_build,
        'vcf_path': str(Path(vcf_path).name),
        'pass_filter': sum(1 for v in variants if v['filter'] == 'PASS'),
        'snvs': sum(1 for v in variants if v['var_type'] == 'SNV'),
        'indels': sum(1 for v in variants if v['var_type'] == 'INDEL'),
        'hom': sum(1 for v in variants if v['zygosity'] == 'HOM'),
        'het': sum(1 for v in variants if v['zygosity'] == 'HET'),
        'high_impact': sum(1 for v in variants if v['impact'] == 'HIGH'),
        'moderate_impact': sum(1 for v in variants if v['impact'] == 'MODERATE'),
        'pathogenic': sum(1 for v in variants if v['classification'] == 'Pathogenic'),
        'likely_pathogenic': sum(1 for v in variants if v['classification'] == 'Likely Pathogenic'),
        'vus': sum(1 for v in variants if v['classification'] == 'VUS'),
        'likely_benign': sum(1 for v in variants if v['classification'] == 'Likely Benign'),
        'benign': sum(1 for v in variants if v['classification'] == 'Benign'),
    }

    # Ti/Tv
    ti_pairs = {frozenset('AG'), frozenset('CT')}
    ti = sum(1 for v in variants if v['var_type'] == 'SNV' and frozenset([v['ref'].upper(), v['alt'].upper()]) in ti_pairs)
    tv = qc['snvs'] - ti
    qc['titv'] = round(ti / tv, 3) if tv > 0 else None

    avg_dp = [int(v['dp']) for v in variants if v['dp'] != '.']
    avg_gq = [int(v['gq']) for v in variants if v['gq'] != '.']
    qc['mean_dp'] = round(sum(avg_dp) / len(avg_dp), 1) if avg_dp else None
    qc['mean_gq'] = round(sum(avg_gq) / len(avg_gq), 1) if avg_gq else None
    qc['het_hom_ratio'] = round(qc['het'] / qc['hom'], 2) if qc['hom'] > 0 else None

    # Genes by variant count
    gene_counts = defaultdict(int)
    for v in variants:
        if v['gene']:
            gene_counts[v['gene']] += 1
    qc['top_genes'] = sorted(gene_counts.items(), key=lambda x: -x[1])[:20]

    return sample_name, variant_map, qc


# ─── PHARMCAT JSON PARSING ───────────────────────────────────────────────────

PHENO_COLORS = {
    'Poor Metabolizer': ('#fef2f2', '#ef4444', '#b91c1c'),
    'Intermediate Metabolizer': ('#fffbeb', '#f59e0b', '#92400e'),
    'Normal Metabolizer': ('#f0fdf4', '#22c55e', '#15803d'),
    'Rapid Metabolizer': ('#eff6ff', '#3b82f6', '#1d4ed8'),
    'Ultrarapid Metabolizer': ('#faf5ff', '#a855f7', '#7e22ce'),
    'Decreased Function': ('#fffbeb', '#f59e0b', '#92400e'),
    'Normal Function': ('#f0fdf4', '#22c55e', '#15803d'),
    'No Function': ('#fef2f2', '#ef4444', '#b91c1c'),
    'No Result': ('#f9fafb', '#9ca3af', '#4b5563'),
    'Indeterminate': ('#faf5ff', '#a855f7', '#7e22ce'),
    'Unknown': ('#f9fafb', '#9ca3af', '#4b5563'),
}

CLASSIFICATION_COLORS = {
    'Strong': ('#fef2f2', '#ef4444'),
    'Moderate': ('#fffbeb', '#f59e0b'),
    'Optional': ('#eff6ff', '#3b82f6'),
    'No recommendation': ('#f9fafb', '#9ca3af'),
    'Unspecified': ('#f9fafb', '#9ca3af'),
}


def pheno_style(phenotype: str) -> tuple[str, str, str]:
    for key, colors in PHENO_COLORS.items():
        if key.lower() in phenotype.lower():
            return colors
    return ('#f9fafb', '#9ca3af', '#4b5563')


def parse_pharmcat_report(report_path: str) -> dict:
    with open(report_path) as f:
        report = json.load(f)

    sample_id = report.get('title', 'Unknown')
    pharmcat_version = report.get('pharmcatVersion', '')
    data_version = report.get('dataVersion', '')
    timestamp = report.get('timestamp', '')

    # ── Genes ──
    genes = {}
    for gene_sym, gdata in report.get('genes', {}).items():
        dips = gdata.get('sourceDiplotypes', [])
        rec_dips = gdata.get('recommendationDiplotypes', dips)

        diplotype_label = 'Unknown'
        phenotypes = []
        activity_score = None
        allele1_fn = allele2_fn = ''
        allele1_name = allele2_name = ''

        if dips:
            d = dips[0]
            diplotype_label = d.get('label', 'Unknown')
            phenotypes = d.get('phenotypes', [])
            activity_score = d.get('activityScore')
            a1 = d.get('allele1') or {}
            a2 = d.get('allele2') or {}
            allele1_name = a1.get('name', '') if isinstance(a1, dict) else ''
            allele2_name = a2.get('name', '') if isinstance(a2, dict) else ''
            allele1_fn = a1.get('function', '') if isinstance(a1, dict) else ''
            allele2_fn = a2.get('function', '') if isinstance(a2, dict) else ''

        pheno_str = ', '.join(phenotypes) if phenotypes else 'No Result'
        bg, fg, border = pheno_style(pheno_str)

        messages = [m.get('message', '') for m in gdata.get('messages', []) if m.get('message')]
        related_drugs = [d.get('name', '') for d in gdata.get('relatedDrugs', [])]
        uncalled = gdata.get('uncalledHaplotypes', [])

        # Calculate coverage from variants list
        vars_list = gdata.get('variants', [])
        total_count = len(vars_list)
        called_count = sum(1 for v in vars_list if v.get('call') and v.get('call') != './.')

        genes[gene_sym] = {
            'gene': gene_sym,
            'chr': gdata.get('chr', ''),
            'diplotype': diplotype_label,
            'phenotypes': phenotypes,
            'pheno_str': pheno_str,
            'pheno_bg': bg, 'pheno_fg': fg, 'pheno_border': border,
            'activity_score': activity_score,
            'allele1_name': allele1_name, 'allele1_fn': allele1_fn,
            'allele2_name': allele2_name, 'allele2_fn': allele2_fn,
            'call_source': gdata.get('callSource', ''),
            'phased': gdata.get('effectivelyPhased', False),
            'messages': messages,
            'related_drugs': related_drugs,
            'uncalled_haplotypes': uncalled,
            'has_undocumented': gdata.get('hasUndocumentedVariations', False),
            'called': called_count,
            'total': total_count
        }

    # Also include unannotated gene calls
    for gdata in report.get('unannotatedGeneCalls', []):
        gene_sym = gdata.get('geneSymbol', '')
        if gene_sym and gene_sym not in genes:
            dips = gdata.get('sourceDiplotypes', [])
            diplotype_label = dips[0].get('label', 'Unknown') if dips else 'Unknown'
            # Unannotated calls often have fewer variant details, but we can still count them
            uvars = gdata.get('variants', [])
            utotal = len(uvars)
            ucalled = sum(1 for v in uvars if v.get('call') and v.get('call') != './.')

            genes[gene_sym] = {
                'gene': gene_sym,
                'chr': gdata.get('chr', ''),
                'diplotype': diplotype_label,
                'phenotypes': [],
                'pheno_str': 'No Recommendation',
                'pheno_bg': '#f9fafb', 'pheno_fg': '#9ca3af', 'pheno_border': '#4b5563',
                'activity_score': None,
                'allele1_name': '', 'allele1_fn': '',
                'allele2_name': '', 'allele2_fn': '',
                'call_source': gdata.get('callSource', ''),
                'phased': gdata.get('effectivelyPhased', False),
                'messages': [],
                'related_drugs': [],
                'uncalled_haplotypes': [],
                'has_undocumented': gdata.get('hasUndocumentedVariations', False),
                'called': ucalled,
                'total': utotal
            }

    # ── Drug Recommendations ──
    drugs = []
    seen = set()
    source_priority = {'CPIC_GUIDELINE': 0, 'DPWG_GUIDELINE': 1, 'FDA_LABEL': 2, 'FDA_ASSOC': 3}

    all_drug_entries = []
    for cat_label, cat_drugs in report.get('drugs', {}).items():
        for drug_name, drug_data in cat_drugs.items():
            source = drug_data.get('source', '')
            for guideline in drug_data.get('guidelines', []):
                for ann in guideline.get('annotations', []):
                    rec = (ann.get('drugRecommendation') or '').strip()
                    if not rec:
                        continue
                    impl = ann.get('implications', [])
                    classification = ann.get('classification', 'Unspecified')
                    phenotypes_used = ann.get('phenotypes', {})
                    dosing = ann.get('dosingInformation', False)
                    alt_drug = ann.get('alternateDrugAvailable', False)
                    other_guidance = ann.get('otherPrescribingGuidance', False)
                    population = ann.get('population', 'general')
                    urls = drug_data.get('urls', [])
                    citations = drug_data.get('citations', [])

                    genes_involved = list(phenotypes_used.keys()) if isinstance(phenotypes_used, dict) else []

                    key = (drug_name, source, classification)
                    priority = source_priority.get(source, 9)

                    all_drug_entries.append({
                        'name': drug_name,
                        'source': source,
                        'source_label': cat_label,
                        'classification': classification,
                        'recommendation': rec,
                        'implications': impl,
                        'phenotypes': phenotypes_used,
                        'genes_involved': genes_involved,
                        'dosing_info': dosing,
                        'alt_drug': alt_drug,
                        'other_guidance': other_guidance,
                        'population': population,
                        'urls': urls,
                        'citations': [c for c in citations if c.get('pmid')],
                        'priority': priority,
                        'key': key,
                    })

    # Deduplicate: keep highest-priority source per drug+classification
    seen_keys = {}
    for entry in sorted(all_drug_entries, key=lambda x: x['priority']):
        dk = (entry['name'], entry['classification'])
        if dk not in seen_keys:
            seen_keys[dk] = entry
            drugs.append(entry)

    # Sort: actionable first (dosing > alt drug), then by classification priority
    cls_order = {'Strong': 0, 'Moderate': 1, 'Optional': 2, 'Unspecified': 3, 'No recommendation': 4}
    drugs.sort(key=lambda d: (
        0 if d['dosing_info'] or d['alt_drug'] else 1,
        cls_order.get(d['classification'], 9),
        d['name']
    ))

    # ── Extract Variants ──
    variants = []
    allele_meta = {}
    for gene_sym, gdata in report.get('genes', {}).items():
        # Map allele names to their metadata (function, activity) for this gene
        for d in gdata.get('sourceDiplotypes', []):
            for a_key in ('allele1', 'allele2'):
                a = d.get(a_key)
                if isinstance(a, dict) and a.get('name'):
                    allele_meta[(gene_sym, a['name'])] = {
                        'function': a.get('function', ''),
                        'activity': a.get('activityValue', '')
                    }

        for v in gdata.get('variants', []):
            chrom = v.get('chromosome', '')
            pos = v.get('position')
            if not chrom or not pos:
                continue
            
            v_alleles = v.get('alleles', [])
            
            # Construct "Related Alleles and Function" and collect activities
            ann_parts = []
            activities = []
            for al in v_alleles:
                meta = allele_meta.get((gene_sym, al))
                if meta:
                    ann_parts.append(f"{al} - {meta['function']}")
                    if meta['activity'] and meta['activity'] != 'n/a':
                        activities.append(str(meta['activity']))
                else:
                    ann_parts.append(al)

            variants.append({
                'gene': gene_sym,
                'chr': chrom,
                'pos': pos,
                'rsid': v.get('dbSnpId') or '',
                'related_ann': '; '.join(ann_parts) if ann_parts else '',
                'activity': ', '.join(set(activities)) if activities else '—',
                'phased': v.get('phased', False),
                'ref': v.get('referenceAllele') or '',
                'call': v.get('call') or '',
            })

    return {
        'sample_id': sample_id,
        'pharmcat_version': pharmcat_version,
        'data_version': data_version,
        'timestamp': timestamp,
        'genes': genes,
        'drugs': drugs,
        'variants': variants,
    }


# ─── HTML REPORT GENERATION ─────────────────────────────────────────────────

def classification_badge_html(cls: str) -> str:
    colors = {
        'Pathogenic':        ('bg-red-100',    'text-red-700',    'border-red-300'),
        'Likely Pathogenic': ('bg-amber-100',  'text-amber-700',  'border-amber-300'),
        'VUS':               ('bg-purple-100', 'text-purple-700', 'border-purple-300'),
        'Likely Benign':     ('bg-green-100',  'text-green-700',  'border-green-300'),
        'Benign':            ('bg-green-50',   'text-green-600',  'border-green-200'),
    }
    bg, text, border = colors.get(cls, ('bg-gray-100', 'text-gray-600', 'border-gray-300'))
    return f'<span class="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium border {bg} {text} {border}">{cls}</span>'


def impact_badge_html(impact: str) -> str:
    colors = {
        'HIGH':     'bg-red-100 text-red-700',
        'MODERATE': 'bg-amber-100 text-amber-700',
        'LOW':      'bg-blue-100 text-blue-700',
        'MODIFIER': 'bg-gray-100 text-gray-500',
    }
    cls = colors.get(impact, 'bg-gray-100 text-gray-500')
    return f'<span class="inline-flex items-center px-1.5 py-0.5 rounded text-xs {cls}">{impact}</span>'


def consequence_display(csq: str) -> str:
    return csq.replace('_', ' ').replace('&', ' + ')


def af_display_html(af: float) -> str:
    label = f'<0.01%' if af < 0.0001 else f'{af*100:.2f}%'
    return f'<span class="text-xs text-gray-600 font-mono whitespace-nowrap">{label}</span>'


def zygosity_badge(zyg: str) -> str:
    if zyg == 'HOM':
        return '<span class="text-amber-600 font-semibold text-xs">HOM</span>'
    return '<span class="text-blue-600 text-xs">HET</span>'


def classification_label(c: str) -> str:
    cls_order = {'Strong': 0, 'Moderate': 1, 'Optional': 2, 'Unspecified': 3, 'No recommendation': 4}
    colors = {
        'Strong':            'bg-red-100 text-red-700 border-red-200',
        'Moderate':          'bg-amber-100 text-amber-700 border-amber-200',
        'Optional':          'bg-blue-100 text-blue-700 border-blue-200',
        'No recommendation': 'bg-gray-100 text-gray-500 border-gray-200',
        'Unspecified':       'bg-gray-100 text-gray-500 border-gray-200',
    }
    cls = colors.get(c, 'bg-gray-100 text-gray-500 border-gray-200')
    return f'<span class="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium border {cls}">{c}</span>'


def source_badge(source: str, label: str) -> str:
    colors = {
        'CPIC_GUIDELINE': 'bg-indigo-50 text-indigo-700',
        'DPWG_GUIDELINE': 'bg-teal-50 text-teal-700',
        'FDA_LABEL':      'bg-orange-50 text-orange-700',
        'FDA_ASSOC':      'bg-pink-50 text-pink-700',
    }
    cls = colors.get(source, 'bg-gray-50 text-gray-600')
    short = {'CPIC_GUIDELINE': 'CPIC', 'DPWG_GUIDELINE': 'DPWG',
             'FDA_LABEL': 'FDA Label', 'FDA_ASSOC': 'FDA Assoc'}
    return f'<span class="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium {cls}">{short.get(source, source)}</span>'


def render_variants_table(variants: list[dict]) -> str:
    """Render a filterable HTML table of all variants."""
    rows = []
    for i, v in enumerate(variants):
        gene_cell = f"""<button onclick="goToGene('{v['gene']}')" class="italic text-indigo-700 font-medium hover:underline text-left">{v['gene']}</button>""" if v['gene'] else '<span class="text-gray-300">—</span>'
        zyg = zygosity_badge(v['zygosity'])
        ref_disp = (v['ref'][:7] + '…') if len(v['ref']) > 10 else v['ref']
        
        # Highlight non-ref alleles
        # If call is different from ref (e.g. call=G/A, ref=G)
        is_variant = (v['zygosity'] != 'REF')
        
        cell_class = "bg-amber-50 font-bold text-amber-900" if is_variant else ""
        
        # Phasing badge
        phasing_html = '<span class="text-indigo-600" title="Phased">🔗</span>' if v.get('phased') else '<span class="text-gray-300">—</span>'
        
        call_display = v['call'] if v['call'] and v['call'] != './.' else '<span class="text-gray-400 italic">No call</span>'
        
        rows.append(
            f'<tr class="hover:bg-gray-50 border-b border-gray-100 variant-row text-xs"'
            f' data-gene="{v["gene"]}" data-zyg="{v["zygosity"]}">'
            f'<td class="px-2 py-2 text-gray-500 font-mono whitespace-nowrap">{v["chr"]}:{v["pos"]}</td>'
            f'<td class="px-2 py-2 font-mono text-xs text-indigo-600">{v["rsid"]}</td>'
            f'<td class="px-2 py-2 font-mono whitespace-nowrap {cell_class}">{call_display}</td>'
            f'<td class="px-2 py-2 font-mono whitespace-nowrap text-center">{phasing_html}</td>'
            f'<td class="px-2 py-2 font-mono whitespace-nowrap">{v["ref"]}</td>'
            f'<td class="px-2 py-2 font-medium whitespace-nowrap">{gene_cell}</td>'
            f'<td class="px-2 py-2 text-[10px] text-gray-700 leading-tight" title="{v["related_ann"]}">{v["related_ann"]}</td>'
            f'<td class="px-2 py-2 whitespace-nowrap">{zyg}</td>'
            f'</tr>'
        )
    return '\n'.join(rows)


def render_gene_cards(pgx_data: dict) -> str:
    if not pgx_data:
        return '<p class="text-gray-400 text-sm">No PharmCAT gene data available.</p>'

    cards = []
    for gene_sym, g in sorted(pgx_data.items()):
        msgs_html = ''
        if g['messages']:
            msgs_html = ''.join(
                f'<p class="text-xs text-amber-700 bg-amber-50 rounded px-2 py-1 mt-1">'
                f'⚠ {m}</p>' for m in g['messages']
            )

        drugs_html = ''
        if g['related_drugs']:
            pills = ' '.join(
                f'<span class="inline-block bg-gray-100 text-gray-600 text-xs px-2 py-0.5 rounded">{d}</span>'
                for d in g['related_drugs']
            )
            drugs_html = f'<div class="mt-2 flex flex-wrap gap-1">{pills}</div>'

        phased_icon = '🔗' if g['phased'] else '❓'
        undoc_warn = '<span class="text-xs text-red-500 ml-1">⚡ undocumented variants</span>' if g['has_undocumented'] else ''

        act_score = ''
        if g['activity_score'] is not None:
            act_score = f'<div class="mt-1 text-xs text-gray-500">Activity score: <span class="font-mono font-medium">{g["activity_score"]}</span></div>'

        cards.append(
            f'<div class="bg-white border border-gray-200 rounded-xl overflow-hidden shadow-sm gene-card" data-gene="{gene_sym}">'
            f'  <div class="bg-gray-50 border-b border-gray-100 px-4 py-3 flex items-start justify-between gap-2">'
            f'    <div>'
            f"""      <div class="font-semibold text-gray-900 text-base italic"><button onclick="goToGene('{gene_sym}')" class="hover:underline text-left">{gene_sym}</button></div>"""
            f'      <div class="text-xs text-gray-500 mt-0.5 font-mono">{g["diplotype"]}</div>'
            f'    </div>'
            f'    <span class="inline-flex items-center px-2.5 py-1 rounded-full text-xs font-medium border mt-0.5"'
            f'     style="background:{g["pheno_bg"]};color:{g["pheno_fg"]};border-color:{g["pheno_border"]}">'
            f'      {g["pheno_str"]}'
            f'    </span>'
            f'  </div>'
            f'  <div class="px-4 py-3">'
            f'    <div class="grid grid-cols-2 gap-x-4 text-xs mb-2">'
            f'      <div class="text-gray-500">Allele 1</div><div class="text-gray-500">Allele 2</div>'
            f'      <div class="font-mono font-medium text-gray-800">{g["allele1_name"] or "—"}</div>'
            f'      <div class="font-mono font-medium text-gray-800">{g["allele2_name"] or "—"}</div>'
            f'      <div class="text-gray-400">{g["allele1_fn"] or ""}</div>'
            f'      <div class="text-gray-400">{g["allele2_fn"] or ""}</div>'
            f'    </div>'
            f'    {act_score}'
            f'    <div class="text-xs text-gray-400 mt-2">{phased_icon} {g["call_source"]} · {g["chr"]}{undoc_warn}</div>'
            f'    {msgs_html}'
            f'    {drugs_html}'
            f'  </div>'
            f'</div>'
        )
    return '\n'.join(cards)


def render_drug_cards(drugs: list[dict]) -> str:
    if not drugs:
        return '<p class="text-gray-400 text-sm">No drug recommendations available.</p>'

    cards = []
    for d in drugs:
        cls_lbl = classification_label(d['classification'])
        src_lbl = source_badge(d['source'], d['source_label'])

        flags = []
        if d['dosing_info']:
            flags.append('<span class="text-xs bg-red-50 text-red-600 border border-red-200 rounded px-1.5 py-0.5">💊 Dosing guidance</span>')
        if d['alt_drug']:
            flags.append('<span class="text-xs bg-orange-50 text-orange-600 border border-orange-200 rounded px-1.5 py-0.5">⇄ Alt drug available</span>')
        if d['other_guidance']:
            flags.append('<span class="text-xs bg-blue-50 text-blue-600 border border-blue-200 rounded px-1.5 py-0.5">ℹ Prescribing guidance</span>')
        flags_html = ' '.join(flags)

        impl_html = ''
        if d['implications']:
            impl_items = ''.join(f'<li>{i}</li>' for i in d['implications'])
            impl_html = f'<ul class="text-xs text-gray-600 mt-2 space-y-0.5 list-disc list-inside">{impl_items}</ul>'

        pheno_html = ''
        if isinstance(d['phenotypes'], dict) and d['phenotypes']:
            pills = ' '.join(
                f"""<button onclick="goToGene('{gene}')" class="font-mono text-xs bg-indigo-50 text-indigo-700 rounded px-1.5 py-0.5 hover:bg-indigo-100 transition-colors">"""
                f"<i>{gene}</i>: {pheno}</button>"
                for gene, pheno in d['phenotypes'].items()
            )
            pheno_html = f'<div class="flex flex-wrap gap-1 mt-2">{pills}</div>'

        citations_html = ''
        if d['citations']:
            cite_links = ' '.join(
                f'<a href="https://pubmed.ncbi.nlm.nih.gov/{c["pmid"]}" target="_blank" '
                f'class="text-xs text-indigo-500 hover:underline">PMID:{c["pmid"]}</a>'
                for c in d['citations'][:3]
            )
            citations_html = f'<div class="mt-2 flex flex-wrap gap-2">{cite_links}</div>'

        url_html = ''
        if d['urls']:
            url_html = (
                f'<a href="{d["urls"][0]}" target="_blank" '
                f'class="text-xs text-indigo-400 hover:underline">Full guideline ↗</a>'
            )

        # Highlight border for actionable drugs
        border_cls = 'border-red-200' if d['dosing_info'] or d['alt_drug'] else 'border-gray-200'
        header_cls = 'bg-red-50' if d['dosing_info'] or d['alt_drug'] else 'bg-gray-50'

        cards.append(
            f'<div class="bg-white border {border_cls} rounded-xl overflow-hidden shadow-sm drug-card" data-genes="{",".join(d["genes_involved"])}">'
            f'  <div class="{header_cls} border-b border-gray-100 px-4 py-3 flex items-start justify-between gap-2">'
            f'    <div>'
            f'      <div class="font-semibold text-gray-900 capitalize">{d["name"]}</div>'
            f'      <div class="flex items-center gap-1.5 mt-1 flex-wrap">'
            f'        {cls_lbl} {src_lbl}'
            f'      </div>'
            f'    </div>'
            f'    <div class="flex flex-wrap gap-1 justify-end">{flags_html}</div>'
            f'  </div>'
            f'  <div class="px-4 py-3">'
            f'    {pheno_html}'
            f'    {impl_html}'
            f'    <p class="text-sm text-gray-800 mt-3 leading-relaxed">{d["recommendation"]}</p>'
            f'    <div class="mt-3 flex items-center gap-3 flex-wrap">'
            f'      {citations_html} {url_html}'
            f'    </div>'
            f'  </div>'
            f'</div>'
        )
    return '\n'.join(cards)


def render_qc_section(qc: dict) -> str:
    def stat(label, val, note='', color='text-gray-900'):
        return (
            f'<div class="bg-white border border-gray-200 rounded-lg px-4 py-3">'
            f'  <div class="text-2xl font-bold font-mono {color}">{val if val is not None else "—"}</div>'
            f'  <div class="text-xs text-gray-500 mt-0.5">{label}</div>'
            f'  {"<div class=\"text-xs text-gray-400\">" + note + "</div>" if note else ""}'
            f'</div>'
        )

    path_color = 'text-red-600' if qc['pathogenic'] > 0 else 'text-gray-900'
    total_str = f'{qc["total"]:,}'
    snv_str = f'{qc["snvs"]:,}'
    indel_str = f'{qc["indels"]:,}'
    pass_str = f'{qc["pass_filter"]:,}'
    pass_pct = f'{round(qc["pass_filter"]/qc["total"]*100,1)}%' if qc["total"] else ""
    dp_str = f'{qc["mean_dp"]}x' if qc["mean_dp"] else None
    parts = [
        stat("Total Variants", total_str, color="text-indigo-600"),
        stat("Pathogenic", qc["pathogenic"], color=path_color),
        stat("LP + VUS", qc["likely_pathogenic"] + qc["vus"]),
        stat("SNVs", snv_str),
        stat("INDELs", indel_str),
        stat("PASS filter", pass_str, pass_pct),
        stat("Mean Depth", dp_str),
        stat("Mean GQ", qc["mean_gq"]),
        stat("Ti/Tv Ratio", qc["titv"]),
        stat("HET/HOM Ratio", qc["het_hom_ratio"]),
    ]
    return '<div class="grid grid-cols-2 sm:grid-cols-3 lg:grid-cols-5 gap-3">' + ''.join(f'  {p}' for p in parts) + '</div>'




def _pheno_pills_html(phenotypes: dict) -> str:
    if not isinstance(phenotypes, dict) or not phenotypes:
        return ''
    pills = ' '.join(
        f"""<button onclick="goToGene('{gene}')" class="font-mono text-xs bg-indigo-50 text-indigo-700 rounded px-1.5 py-0.5 hover:bg-indigo-100 transition-colors">"""
        f"<i>{gene}</i>: {pheno}</button>"
        for gene, pheno in phenotypes.items()
    )
    return f'<div class="flex flex-wrap gap-1 mb-2">{pills}</div>'


def _render_actionable_drugs(pgx: dict | None) -> str:
    if not pgx:
        return ''
    cards = []
    for d in pgx['drugs']:
        if not (d['dosing_info'] or d['alt_drug']):
            continue
        pheno_html = _pheno_pills_html(d['phenotypes'])
        cls_lbl = classification_label(d['classification'])
        src_lbl = source_badge(d['source'], d['source_label'])

        # Implications
        impl_html = ''
        if d.get('implications'):
            items = ''.join(f'<li>{i}</li>' for i in d['implications'])
            impl_html = f'<ul class="text-xs text-gray-600 mt-2 space-y-0.5 list-disc list-inside">{items}</ul>'

        # Actionable flags
        flags = []
        if d['dosing_info']:
            flags.append('<span class="text-[10px] bg-red-100 text-red-700 px-1.5 py-0.5 rounded border border-red-200">💊 Dosing guidance</span>')
        if d['alt_drug']:
            flags.append('<span class="text-[10px] bg-orange-100 text-orange-700 px-1.5 py-0.5 rounded border border-orange-200">⇄ Alt drug available</span>')
        flags_html = ' '.join(flags)

        url_html = ''
        if d.get('urls'):
            url_html = f'<a href="{d["urls"][0]}" target="_blank" class="text-xs text-indigo-400 hover:underline">Full guideline ↗</a>'

        card = (
            f'<div class="bg-white border border-red-200 rounded-xl overflow-hidden shadow-sm hover:shadow-md transition-shadow drug-card" data-genes="{",".join(d["genes_involved"])}">'
            f'  <div class="bg-red-50 border-b border-red-100 px-4 py-3 flex items-center justify-between gap-2">'
            f'    <div class="flex items-center gap-2">'
            f'      <div class="font-semibold capitalize text-gray-900">{d["name"]}</div>'
            f'      {flags_html}'
            f'    </div>'
            f'    <div class="flex gap-1">{cls_lbl} {src_lbl}</div>'
            f'  </div>'
            f'  <div class="px-4 py-3">'
            f'    {pheno_html}'
            f'    {impl_html}'
            f'    <p class="text-sm text-gray-800 mt-3 leading-relaxed font-medium">{d["recommendation"]}</p>'
            f'    <div class="mt-3">{url_html}</div>'
            f'  </div>'
            f'</div>'
        )
        cards.append(card)
    if not cards:
        return '<p class="text-gray-400 text-sm">No dosing-specific recommendations found.</p>'
    return '\n'.join(cards)

def build_html(vcf_path: str, report_path: str, sample_name: str,
               variants: list[dict], qc: dict, pgx: dict | None) -> str:

    generated_at = datetime.now().strftime('%Y-%m-%d %H:%M')
    genes_html = render_gene_cards(pgx['genes'] if pgx else {})
    drugs_html = render_drug_cards(pgx['drugs'] if pgx else [])
    MAX_TABLE_ROWS = 5000
    if len(variants) > MAX_TABLE_ROWS:
        variants_truncated = sorted(
            variants,
            key=lambda v: ({'HIGH':0,'MODERATE':1,'LOW':2,'MODIFIER':3}.get(v['impact'],9),
                           {'Pathogenic':0,'Likely Pathogenic':1,'VUS':2,'Likely Benign':3,'Benign':4}.get(v['classification'],5))
        )[:MAX_TABLE_ROWS]
        truncation_note = f'<tr><td colspan="12" class="text-center text-amber-700 bg-amber-50 py-3 text-xs">⚠ Showing top {MAX_TABLE_ROWS:,} variants by impact/classification. {len(variants):,} total in VCF. Export TSV for full list.</td></tr>'
    else:
        variants_truncated = variants
        truncation_note = ''
    variants_html = render_variants_table(variants_truncated) + truncation_note
    qc_html = render_qc_section(qc)

    pgx_sample_id = pgx['sample_id'] if pgx else sample_name
    pharmcat_version = pgx['pharmcat_version'] if pgx else '—'
    data_version = pgx['data_version'] if pgx else '—'

    n_genes = len(pgx['genes']) if pgx else 0
    n_drugs = len(pgx['drugs']) if pgx else 0
    n_actionable = sum(1 for d in (pgx['drugs'] if pgx else []) if d['dosing_info'] or d['alt_drug'])
    n_nonref = sum(1 for v in variants if v.get('zygosity') != 'REF')

    # High/moderate impact variants for summary table
    # Significant Variants for summary table
    # Instead of just Pathogenic (which is rare in PGx), show any variant with a non-reference call
    # that has a PharmCAT related annotation.
    def is_variant_call(v):
        call = v['call'] or './.'
        if call == './.' or not v['ref']: return False
        parts = call.replace('|', '/').split('/')
        return any(p != v['ref'] and p != '.' for p in parts)

    sig_variants = [v for v in variants if is_variant_call(v) and v.get('related_ann')]
    sig_rows = render_variants_table(sig_variants) if sig_variants else '<tr><td colspan="12" class="text-center text-gray-400 py-6 text-sm">No non-reference clinical variants detected</td></tr>'

    return f'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>NXF-ALIGNMENT — {pgx_sample_id}</title>
<script src="https://cdn.tailwindcss.com"></script>
<style>
  @import url('https://fonts.googleapis.com/css2?family=IBM+Plex+Mono:wght@400;500&family=Inter:wght@300;400;500;600;700&display=swap');
  body {{ font-family: 'Inter', sans-serif; }}
  .font-mono {{ font-family: 'IBM Plex Mono', monospace; }}
  @media print {{
    .no-print {{ display: none !important; }}
    .page-break {{ page-break-before: always; }}
    body {{ font-size: 11px; }}
  }}
  .tab-content {{ display: none; }}
  .tab-content.active {{ display: block; }}
  .tab-btn.active {{ background-color: #eef2ff; color: #4338ca; font-weight: 600; border-bottom: 2px solid #4338ca; }}
</style>
</head>
<body class="bg-gray-50 text-gray-900 min-h-screen">


<!-- REPORT HEADER BANNER -->
<div class="bg-indigo-700 text-white">
  <div class="max-w-screen-xl mx-auto px-6 py-6">
    <div class="flex flex-wrap items-start justify-between gap-4">
      <div>
        <div class="text-xs font-medium text-indigo-300 uppercase tracking-widest mb-1">Clinical Pharmacogenomics Report</div>
        <div class="text-3xl font-bold">{pgx_sample_id}</div>
        <div class="text-indigo-200 text-sm mt-1">Generated by the NXF-ALIGNMENT pipeline · {generated_at} · Genome build: {qc["genome_build"]}</div>
      </div>
      <div class="grid grid-cols-3 gap-3 text-center">
        <div class="bg-white/10 rounded-lg px-4 py-3">
          <div class="text-2xl font-bold">{n_genes}</div>
          <div class="text-xs text-indigo-200">PGx Genes</div>
        </div>
        <div class="bg-white/10 rounded-lg px-4 py-3">
          <div class="text-2xl font-bold">{n_drugs}</div>
          <div class="text-xs text-indigo-200">Drug Recs</div>
        </div>
        <div class="bg-white/10 rounded-lg px-4 py-3">
          <div class="text-2xl font-bold">{n_actionable}</div>
          <div class="text-xs text-indigo-200">Actionable</div>
        </div>
      </div>
    </div>

    <!-- Meta row -->
    <div class="mt-4 flex flex-wrap gap-4 text-xs text-indigo-300">
      <span>📄 VCF: <span class="text-white font-mono">{qc["vcf_path"]}</span></span>
      <span>🔬 PharmCAT: <span class="text-white">{pharmcat_version}</span></span>
      <span>📚 PharmCAT data: <span class="text-white">{data_version}</span></span>
      <span>🧬 Sample: <span class="text-white font-mono">{sample_name}</span></span>
      <span>📊 Variants: <span class="text-white font-mono">{qc["total"]:,}</span></span>
    </div>
  </div>
</div>

<!-- DISCLAIMER -->
<div class="bg-amber-50 border-b border-amber-200">
  <div class="max-w-screen-xl mx-auto px-6 py-2 text-xs text-amber-800">
    ⚠ <strong>Research use only.</strong> This report is not a substitute for clinical laboratory testing or physician interpretation. 
    Variant classifications are computational predictions and have not been clinically validated.
  </div>
</div>

<!-- NAVIGATION TABS -->
<div class="bg-white border-b border-gray-200 no-print sticky top-0 z-30">
  <div class="max-w-screen-xl mx-auto px-6 flex gap-0">
    <button class="tab-btn active px-4 py-3 text-sm transition-colors" onclick="showTab('summary')">Summary</button>
    <button class="tab-btn px-4 py-3 text-sm text-gray-600 hover:text-gray-900" onclick="showTab('pgx-genes')">PGx Genes ({n_genes})</button>
    <button class="tab-btn px-4 py-3 text-sm text-gray-600 hover:text-gray-900" onclick="showTab('drugs')">Drug Recommendations ({n_drugs})</button>
    <button class="tab-btn px-4 py-3 text-sm text-gray-600 hover:text-gray-900" onclick="showTab('variants')">All Variants ({len(variants):,} total, {n_nonref:,} non-ref)</button>
    <button class="tab-btn px-4 py-3 text-sm text-gray-600 hover:text-gray-900" onclick="showTab(\'qc\')">QC Metrics</button>
  </div>
</div>

<div class="active-filter-info hidden bg-indigo-600 text-white shadow-lg sticky top-[53px] z-20">
  <div class="max-w-screen-xl mx-auto px-6 py-2.5 flex items-center justify-between gap-4">
    <div class="flex items-center gap-3">
      <span class="bg-white/20 p-1.5 rounded-lg text-lg">🧬</span>
      <div>
        <div class="text-[10px] text-indigo-200 uppercase font-bold tracking-wider leading-none mb-0.5">Active Filter</div>
        <div class="text-sm font-medium">Showing information for <span class="active-gene-name font-bold italic underline decoration-indigo-300"></span> gene</div>
      </div>
    </div>
    <button onclick="resetFilters()" class="bg-white/10 hover:bg-white/20 text-white text-xs px-3 py-1.5 rounded-lg border border-white/20 transition-all font-semibold backdrop-blur-sm">Clear Filter</button>
  </div>
</div>

<div class="max-w-screen-xl mx-auto px-6 py-6">

<!-- ════ SUMMARY TAB ════ -->
<div id="tab-summary" class="tab-content active">

  <!-- PGx Gene Summary Table -->
  <section class="mb-6">
    <div class="flex items-center gap-2 mb-3 cursor-pointer group" onclick="toggleSection('summary-genes')">
      <h2 class="text-sm font-semibold text-gray-500 uppercase tracking-widest group-hover:text-indigo-600 transition-colors">Pharmacogenomics — Gene Summary</h2>
      <span id="summary-genes-icon" class="text-gray-400 transition-transform duration-200">▼</span>
    </div>
    <div id="summary-genes" class="bg-white border border-gray-200 rounded-xl overflow-hidden shadow-sm">
      <table class="w-full text-sm">
        <thead>
          <tr class="bg-gray-50 border-b border-gray-200 text-xs text-gray-500 uppercase tracking-wide">
            <th class="px-4 py-3 text-left">Gene</th>
            <th class="px-4 py-3 text-left">Genotypes</th>
            <th class="px-4 py-3 text-left">Diplotype</th>
            <th class="px-4 py-3 text-left">Phenotype</th>
            <th class="px-4 py-3 text-left">Related Drugs</th>
            <th class="px-4 py-3 text-left">Source</th>
          </tr>
        </thead>
        <tbody>
          {''.join(
            f'<tr class="border-b border-gray-100 hover:bg-gray-50 summary-gene-row" data-gene="{g["gene"]}">'
            f"""<td class="px-4 py-2 font-medium italic text-indigo-700"><button onclick="goToGene('{g['gene']}')" class="hover:underline text-left">{g['gene']}</button></td>"""
            f'<td class="px-4 py-2 font-mono text-xs text-gray-500">{g["called"]} / {g["total"]}</td>'
            f'<td class="px-4 py-2 font-mono text-xs">{g["diplotype"]}</td>'
            f'<td class="px-4 py-2">'
            f'  <span class="inline-flex items-center px-2 py-0.5 rounded text-xs font-medium border"'
            f'   style="background:{g["pheno_bg"]};color:{g["pheno_fg"]};border-color:{g["pheno_border"]}">'
            f'    {g["pheno_str"]}</span></td>'
            f'<td class="px-4 py-2 text-xs text-gray-500">{", ".join(g["related_drugs"][:4]) or "—"}</td>'
            f'<td class="px-4 py-2 text-xs text-gray-400">{g["call_source"]}</td>'
            f'</tr>'
            for g in sorted(pgx['genes'].values(), key=lambda x: x['gene'])
          ) or "<tr><td colspan='6' class='text-center text-gray-400 py-6 text-sm'>No PharmCAT gene data</td></tr>"}
        </tbody>
      </table>
    </div>
  </section>

  <!-- Actionable Drug Recommendations -->
  <section class="mb-6">
    <div class="flex items-center gap-2 mb-3 cursor-pointer group" onclick="toggleSection('summary-drugs')">
      <h2 class="text-sm font-semibold text-gray-500 uppercase tracking-widest group-hover:text-indigo-600 transition-colors">
        Actionable Drug Recommendations
        <span class="ml-2 text-xs font-normal bg-red-100 text-red-600 rounded px-1.5 py-0.5 normal-case tracking-normal">
          {n_actionable} with dosing or alt drug guidance
        </span>
      </h2>
      <span id="summary-drugs-icon" class="text-gray-400 transition-transform duration-200">▼</span>
    </div>
    <div id="summary-drugs" class="grid grid-cols-1 lg:grid-cols-2 gap-4">
      {_render_actionable_drugs(pgx)}
    </div>
  </section>

  <!-- Detected Variants with Clinical Annotations -->
  <section class="mb-6">
    <div class="flex items-center gap-2 mb-3 cursor-pointer group" onclick="toggleSection('summary-variants')">
      <h2 class="text-sm font-semibold text-gray-500 uppercase tracking-widest group-hover:text-indigo-600 transition-colors">
        Detected Variants with Clinical Annotations
      </h2>
      <span id="summary-variants-icon" class="text-gray-400 transition-transform duration-200">▼</span>
    </div>
    <div id="summary-variants" class="bg-white border border-gray-200 rounded-xl overflow-hidden shadow-sm">
      <div class="overflow-x-auto">
        <table class="w-full text-sm">
          <thead>
            <tr class="bg-gray-50 border-b border-gray-200 text-[10px] text-gray-500 uppercase tracking-wider">
              <th class="px-2 py-2 text-left font-semibold">Position in VCF</th>
              <th class="px-2 py-2 text-left font-semibold">RSID</th>
              <th class="px-2 py-2 text-left font-semibold">Call in VCF</th>
              <th class="px-2 py-2 text-center font-semibold">Phase</th>
              <th class="px-2 py-2 text-left font-semibold">Reference</th>
              <th class="px-2 py-2 text-left font-semibold">Gene</th>
              <th class="px-2 py-2 text-left font-semibold">Related Alleles and Function</th>
              <th class="px-2 py-2 text-left font-semibold">Zyg</th>
            </tr>
          </thead>
          <tbody>{sig_rows}</tbody>
        </table>
      </div>
    </div>
  </section>


</div><!-- /summary -->

<!-- ════ PGx GENES TAB ════ -->
<div id="tab-pgx-genes" class="tab-content">
  <h2 class="text-sm font-semibold text-gray-500 uppercase tracking-widest mb-4">Pharmacogenomics Gene Calls</h2>
  <div class="grid grid-cols-1 md:grid-cols-2 xl:grid-cols-3 gap-4">
    {genes_html}
  </div>
</div>

<!-- ════ DRUGS TAB ════ -->
<div id="tab-drugs" class="tab-content">
  <div class="flex items-center gap-4 mb-4 no-print">
    <h2 class="text-sm font-semibold text-gray-500 uppercase tracking-widest">Drug Recommendations</h2>
    <div class="flex gap-2 ml-auto">
      <label class="text-xs text-gray-500 flex items-center gap-1">
        <input type="checkbox" id="filterActionable" onchange="filterDrugs()" class="accent-indigo-600">
        Actionable only
      </label>
    </div>
  </div>
  <div id="drugsContainer" class="grid grid-cols-1 lg:grid-cols-2 gap-4">
    {drugs_html}
  </div>
</div>

<!-- ════ VARIANTS TAB ════ -->
<div id="tab-variants" class="tab-content">
  <!-- Filters -->
  <div class="bg-white border border-gray-200 rounded-xl p-4 mb-4 no-print">
    <div class="flex flex-wrap gap-3 items-end">
      <div>
        <label class="block text-xs text-gray-500 mb-1">Gene</label>
        <select id="filterGene" class="border border-gray-200 rounded-lg px-3 py-1.5 text-sm focus:outline-none focus:border-indigo-400" onchange="filterVariants()">
          <option value="">All Genes</option>
          {''.join(f'<option value="{g}">{g}</option>' for g in sorted(list(set(v['gene'] for v in variants if v.get('gene')))))}
        </select>
      </div>
      <div>
        <label class="block text-xs text-gray-500 mb-1">Zygosity</label>
        <select id="filterZyg" class="border border-gray-200 rounded-lg px-3 py-1.5 text-sm focus:outline-none focus:border-indigo-400" onchange="filterVariants()">
          <option value="">All</option>
          <option>HOM</option>
          <option>HET</option>
        </select>
      </div>
      <button onclick="resetFilters()" class="text-xs text-gray-400 hover:text-gray-700 px-3 py-1.5 border border-gray-200 rounded-lg">Reset</button>
      <span id="variantCount" class="text-xs text-gray-400 ml-auto self-center font-mono">{len(variants):,} total ({n_nonref:,} non-ref)</span>
    </div>
  </div>

  <div class="bg-white border border-gray-200 rounded-xl overflow-hidden shadow-sm">
    <div class="overflow-x-auto">
      <table class="w-full text-sm" id="variantsTable">
        <thead>
          <tr class="bg-gray-50 border-b border-gray-200 text-[10px] text-gray-500 uppercase tracking-wider">
            <th class="px-2 py-2 text-left font-semibold">Position in VCF</th>
            <th class="px-2 py-2 text-left font-semibold">RSID</th>
            <th class="px-2 py-2 text-left font-semibold">Call in VCF</th>
            <th class="px-2 py-2 text-center font-semibold">Phase</th>
            <th class="px-2 py-2 text-left font-semibold">Reference</th>
            <th class="px-2 py-2 text-left font-semibold">Gene</th>
            <th class="px-2 py-2 text-left font-semibold">Related Alleles and Function</th>
            <th class="px-2 py-2 text-left font-semibold">Zyg</th>
          </tr>
        </thead>
        <tbody id="variantsBody">
          {variants_html}
        </tbody>
      </table>
    </div>
  </div>
</div>

<!-- ════ QC TAB ════ -->
<div id="tab-qc" class="tab-content">
  <div class="mb-8">
    <h2 class="text-sm font-semibold text-gray-500 uppercase tracking-widest mb-3">Sequencing Quality Overview</h2>
    {qc_html}
  </div>

  <h2 class="text-sm font-semibold text-gray-500 uppercase tracking-widest mb-4">Detailed Quality Metrics</h2>

  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <div class="bg-white border border-gray-200 rounded-xl p-5 shadow-sm">
      <h3 class="font-semibold text-gray-700 mb-4">Variant Summary</h3>
      {''.join(
        f'<div class="flex justify-between py-2 border-b border-gray-50 text-sm">'
        f'  <span class="text-gray-500">{k}</span>'
        f'  <span class="font-mono font-medium">{v}</span>'
        f'</div>'
        for k, v in [
          ("Total variants", f'{qc["total"]:,}'),
          ("PASS filter", f'{qc["pass_filter"]:,} ({round(qc["pass_filter"]/qc["total"]*100,1) if qc["total"] else 0}%)'),
          ("SNVs", f'{qc["snvs"]:,}'),
          ("INDELs", f'{qc["indels"]:,}'),
          ("Ti/Tv ratio", qc["titv"] or "—"),
          ("Heterozygous", f'{qc["het"]:,}'),
          ("Homozygous", f'{qc["hom"]:,}'),
          ("HET/HOM ratio", qc["het_hom_ratio"] or "—"),
          ("Mean depth (DP)", f'{qc["mean_dp"]}×' if qc["mean_dp"] else "—"),
          ("Mean GQ", qc["mean_gq"] or "—"),
        ]
      )}
    </div>

    <div class="bg-white border border-gray-200 rounded-xl p-5 shadow-sm">
      <h3 class="font-semibold text-gray-700 mb-4">Classification Breakdown</h3>
      {''.join(
        f'<div class="flex justify-between py-2 border-b border-gray-50 text-sm items-center">'
        f'  <span>{classification_badge_html(k)}</span>'
        f'  <span class="font-mono text-gray-700">{v} <span class="text-gray-400 text-xs">({round(v/qc["total"]*100,1) if qc["total"] else 0}%)</span></span>'
        f'</div>'
        for k, v in [
          ("Pathogenic", qc["pathogenic"]),
          ("Likely Pathogenic", qc["likely_pathogenic"]),
          ("VUS", qc["vus"]),
          ("Likely Benign", qc["likely_benign"]),
          ("Benign", qc["benign"]),
        ]
      )}
    </div>

    <div class="bg-white border border-gray-200 rounded-xl p-5 shadow-sm">
      <h3 class="font-semibold text-gray-700 mb-4">Impact Distribution</h3>
      {''.join(
        f'<div class="flex justify-between py-2 border-b border-gray-50 text-sm items-center">'
        f'  <span>{impact_badge_html(k)}</span>'
        f'  <span class="font-mono text-gray-700">{v:,}</span>'
        f'</div>'
        for k, v in [
          ("HIGH", qc["high_impact"]),
          ("MODERATE", qc["moderate_impact"]),
          ("LOW", sum(1 for x in variants if x["impact"] == "LOW")),
          ("MODIFIER", sum(1 for x in variants if x["impact"] == "MODIFIER")),
        ]
      )}
    </div>

    <div class="bg-white border border-gray-200 rounded-xl p-5 shadow-sm">
      <h3 class="font-semibold text-gray-700 mb-4">Top Genes by Variant Count</h3>
      {''.join(
        f'<div class="flex justify-between py-1.5 border-b border-gray-50 text-sm">'
        f'  <span class="italic text-indigo-700">{gene}</span>'
        f'  <span class="font-mono text-gray-600">{cnt}</span>'
        f'</div>'
        for gene, cnt in qc["top_genes"][:15]
      )}
    </div>
  </div>
</div>

</div><!-- /main container -->

<!-- FOOTER -->
<div class="border-t border-gray-200 bg-white mt-8 py-4">
  <div class="max-w-screen-xl mx-auto px-6 text-xs text-gray-400 flex justify-between">
    <span>VarPath · Generated {generated_at}</span>
    <span>PharmCAT {pharmcat_version} · {data_version}</span>
  </div>
</div>

<script>
function showTab(name) {{
  document.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
  document.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
  document.getElementById('tab-' + name).classList.add('active');
  
  // Update button state
  const btns = document.querySelectorAll('.tab-btn');
  btns.forEach(btn => {{
    if (btn.getAttribute('onclick').includes("'" + name + "'")) {{
      btn.classList.add('active');
    }}
  }});
}}

function goToGene(gene) {{
  const filterEl = document.getElementById('filterGene');
  if (filterEl) {{
    filterEl.value = gene;
  }}
  applyFilters();
  showTab('variants');
}}

function toggleSection(id) {{
  const content = document.getElementById(id);
  const icon = document.getElementById(id + '-icon');
  if (content.classList.contains('hidden')) {{
    content.classList.remove('hidden');
    if (icon) icon.style.transform = 'rotate(0deg)';
  }} else {{
    content.classList.add('hidden');
    if (icon) icon.style.transform = 'rotate(-90deg)';
  }}
}}

function applyFilters() {{
  const geneVal = document.getElementById('filterGene').value;
  const gene = geneVal.toLowerCase();
  const zyg = document.getElementById('filterZyg').value;
  const actionableOnly = document.getElementById('filterActionable') ? document.getElementById('filterActionable').checked : false;

  // 1. Update Global Banner
  const infoEls = document.querySelectorAll('.active-filter-info');
  const geneNameEls = document.querySelectorAll('.active-gene-name');
  if (geneVal) {{
    infoEls.forEach(el => el.classList.remove('hidden'));
    geneNameEls.forEach(el => el.textContent = geneVal);
  }} else {{
    infoEls.forEach(el => el.classList.add('hidden'));
  }}

  // 2. Filter Variants Table
  let visibleVariants = 0;
  let visibleNonRef = 0;
  document.querySelectorAll('.variant-row').forEach(row => {{
    const show = (!gene || row.dataset.gene.toLowerCase() === gene) &&
                 (!zyg || row.dataset.zyg === zyg);
    row.style.display = show ? '' : 'none';
    if (show) {{
      visibleVariants++;
      if (row.dataset.zyg !== 'REF') visibleNonRef++;
    }}
  }});
  const countEl = document.getElementById('variantCount');
  if (countEl) countEl.textContent = visibleVariants + ' total (' + visibleNonRef + ' non-ref)';

  // 3. Filter Gene Cards (PGx Genes tab)
  document.querySelectorAll('.gene-card').forEach(card => {{
    const show = !gene || card.dataset.gene.toLowerCase() === gene;
    card.style.display = show ? '' : 'none';
  }});

  // 4. Filter Drug Cards (Drug Recommendations tab & Summary tab)
  document.querySelectorAll('.drug-card').forEach(card => {{
    const geneMatch = !gene || (card.dataset.genes && card.dataset.genes.toLowerCase().split(',').includes(gene));
    const actionMatch = !actionableOnly || card.classList.contains('border-red-200');
    card.style.display = (geneMatch && actionMatch) ? '' : 'none';
  }});

  // 5. Filter Summary Gene Table
  document.querySelectorAll('.summary-gene-row').forEach(row => {{
    const show = !gene || row.dataset.gene.toLowerCase() === gene;
    row.style.display = show ? '' : 'none';
  }});
}}

function filterVariants() {{
  applyFilters();
}}

function filterDrugs() {{
  applyFilters();
}}

function resetFilters() {{
  const gEl = document.getElementById('filterGene');
  const zEl = document.getElementById('filterZyg');
  if (gEl) gEl.value = '';
  if (zEl) zEl.value = '';
  applyFilters();
}}
</script>
</body>
</html>'''


# ─── CLI ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='VarPath: VCF + PharmCAT JSON → HTML Clinical Report',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--vcf',    required=True, help='Annotated VCF file (.vcf, .vcf.gz, .bgz)')
    parser.add_argument('--report', required=True, help='PharmCAT report JSON (e.g. sample_report.json)')

    parser.add_argument('--out',    default=None,  help='Output HTML file (default: <sample>_report.html)')
    args = parser.parse_args()

    # ── Parse Reports ──
    print(f'[1/3] Parsing VCF: {args.vcf}')
    sample_name, vcf_map, qc = parse_vcf(args.vcf)
    
    print(f'[2/3] Parsing PharmCAT report: {args.report}')
    pgx = parse_pharmcat_report(args.report)

    # ── Merge Variants ──
    # The table is now driven by PharmCAT variants, supplemented by VCF metrics
    final_variants = []
    for pv in pgx['variants']:
        key = (pv['chr'], pv['pos'])
        vcf_match = vcf_map.get(key)
        
        if not vcf_match:
            # Check without 'chr' prefix
            alt_key = (pv['chr'][3:] if pv['chr'].startswith('chr') else f"chr{pv['chr']}", pv['pos'])
            vcf_match = vcf_map.get(alt_key)

        if vcf_match:
            # Update with technical metrics from VCF
            pv.update({
                'alt': vcf_match['alt'],
                'ref': vcf_match['ref'],
                'var_type': vcf_match['var_type'],
                'consequence': vcf_match['consequence'],
                'impact': vcf_match['impact'],
                'hgvsc': vcf_match['hgvsc'],
                'hgvsp': vcf_match['hgvsp'],
                'zygosity': vcf_match['zygosity'],
                'gt': vcf_match['gt'],
                'dp': vcf_match['dp'],
                'af': vcf_match['af'],
                'gnomad_af': vcf_match['gnomad_af'],
                'classification': vcf_match['classification'],
            })
            # Prefer VCF genotype if available for the "Call" column
            pv['call'] = vcf_match['gt']
        else:
            # Fallback for variants not in VCF (e.g. wild-type positions)
            pv.update({
                'alt': '',
                'var_type': 'SNV',
                'consequence': 'Reference',
                'impact': 'MODIFIER',
                'hgvsc': '',
                'hgvsp': '',
                'zygosity': gt_zygosity(pv['call'] or './.', pv['ref']),
                'gt': pv['call'] or './.',
                'dp': '.',
                'af': 0.0,
                'gnomad_af': 0.0,
                'classification': 'Benign',
            })
        final_variants.append(pv)
    
    variants = final_variants

    # ── Render ──
    print(f'[3/3] Rendering HTML report…')
    html = build_html(args.vcf, args.report, sample_name, variants, qc, pgx)

    out_path = args.out or f'{pgx["sample_id"]}-pgx-report.html'
    Path(out_path).write_text(html, encoding='utf-8')
    print(f'\n✓ Report written to: {out_path}')
    print(f'  Variants: {qc["total"]:,}  |  Pathogenic: {qc["pathogenic"]}  |  Actionable drugs: {sum(1 for d in pgx["drugs"] if d["dosing_info"] or d["alt_drug"])}')


if __name__ == '__main__':
    main()
