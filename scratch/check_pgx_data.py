import json
import sys

path = '/Users/angeloas/code/nxf-alignment/output/04-pgx/reads3_pharmcat/reads3.report.json'
try:
    with open(path) as f:
        r = json.load(f)
except Exception as e:
    print(f"Error loading JSON: {e}")
    sys.exit(1)

genes = r.get('genes', {})
for g, data in genes.items():
    sd = data.get('sourceDiplotypes', [])
    if sd:
        called = sd[0].get('genotypesCalled')
        total = sd[0].get('callsAtPositions')
        print(f"{g}: {called} / {total}")
    else:
        print(f"{g}: No data")
