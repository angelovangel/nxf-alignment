import json
path = '/Users/angeloas/code/nxf-alignment/output/04-pgx/reads3_pharmcat/reads3.report.json'
with open(path) as f:
    r = json.load(f)
genes = r.get('genes', {})
for g_sym, data in genes.items():
    vars = data.get('variants', [])
    total = len(vars)
    called = sum(1 for v in vars if v.get('call'))
    print(f"{g_sym}: {called} / {total}")
