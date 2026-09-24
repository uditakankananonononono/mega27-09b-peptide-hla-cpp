"""Parse EMBOSS pepinfo (EBI REST) raw outputs into per-property residue fractions.

Reads /tmp/pepinfo_raw/<CPP>.out (written by emboss_pepinfo_screen.py poll),
computes the fraction of residues flagged 1 in each property block, and writes
results/cpp_emboss_pepinfo.json with per-CPP and panel-mean fractions.
"""
import re, json, glob, os

props_all, per_cpp = {}, {}
for f in sorted(glob.glob('/tmp/pepinfo_raw/*.out')):
    name = os.path.basename(f).split('.')[0]
    txt = open(f).read()
    blocks = re.findall(r'Printing out (\w[\w ]*?) residues.*?Position  Residue\s*Result\n(.*?)(?=\n\n\n|\Z)', txt, re.S)
    fracs = {}
    for prop, body in blocks:
        vals = [int(l.split()[-1]) for l in body.strip().splitlines() if l.split() and l.split()[-1] in ('0', '1')]
        if vals:
            fracs[prop.strip()] = round(sum(vals) / len(vals), 3)
    per_cpp[name] = fracs
    for p, v in fracs.items():
        props_all.setdefault(p, []).append(v)
summary = {p: round(sum(v) / len(v), 3) for p, v in sorted(props_all.items())}
json.dump({'n_cpps': len(per_cpp), 'mean_fraction_per_property': summary, 'per_cpp': per_cpp},
          open('results/cpp_emboss_pepinfo.json', 'w'), indent=1)
print(len(per_cpp), summary)
