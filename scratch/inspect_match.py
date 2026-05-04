import json
path = '/Users/angeloas/code/nxf-alignment/output/04-pgx/reads3_pharmcat/reads3.match.json'
with open(path) as f:
    r = json.load(f)
genes = r.get('genes', {})
for g in ['CACNA1S', 'DPYD', 'UGT1A1']:
    print(f"--- {g} ---")
    data = genes.get(g, {})
    # print keys and some interesting values
    print(f"Keys: {list(data.keys())}")
    if 'results' in data:
        # PharmCAT match.json usually has results
        res = data['results']
        if isinstance(res, list) and len(res) > 0:
            match = res[0]
            print(f"genotypesCalled: {match.get('genotypesCalled')}")
            print(f"callsAtPositions: {match.get('callsAtPositions')}")
