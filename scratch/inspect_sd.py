import json
path = '/Users/angeloas/code/nxf-alignment/output/04-pgx/reads3_pharmcat/reads3.report.json'
with open(path) as f:
    r = json.load(f)
genes = r.get('genes', {})
for g in ['CACNA1S', 'DPYD', 'UGT1A1']:
    sd = genes.get(g, {}).get('sourceDiplotypes', [])
    if sd:
        print(f"--- {g} ---")
        print(json.dumps(sd[0], indent=2))
