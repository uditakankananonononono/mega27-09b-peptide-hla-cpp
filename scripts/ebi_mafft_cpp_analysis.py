"""Parse EBI MAFFT (REST) outputs for the 18 designed CPPs + 7 CPPsite archetypes.

Genuine external tool use: EBI MAFFT REST job (multiple sequence alignment,
guide tree, percent identity matrix). Descriptive family-level novelty check:
how much pairwise identity do designed CPPs share with known CPP archetypes?
Ultra-short compositionally-biased peptides => alignment is descriptive only.
"""
import json, re

pim = open('results/ebi_mafft_cpp_pim.txt').read().splitlines()
rows = []
for ln in pim:
    m = re.match(r'\s*\d+:\s+(\S+)\s+(.*)', ln)
    if m:
        rows.append((m.group(1), [float(x) for x in m.group(2).split()]))
names = [r[0] for r in rows]
arch_names = [n for n in names if n.startswith('CPPsite')]
des_names = [n for n in names if n.startswith('9B-CPP')]
mat = {n: dict(zip(names, vals)) for n, vals in rows}

per_designed = {}
for d in des_names:
    vs = {a: mat[d][a] for a in arch_names}
    nearest = max(vs, key=vs.get)
    per_designed[d] = {
        'identity_vs_archetypes': vs,
        'nearest_archetype': nearest,
        'nearest_identity': vs[nearest],
        'mean_identity_vs_archetypes': round(sum(vs.values()) / len(vs), 2),
    }

within = [mat[a][b] for i, a in enumerate(des_names) for b in des_names[i+1:]]
cross = [mat[d][a] for d in des_names for a in arch_names]
arch_within = [mat[a][b] for i, a in enumerate(arch_names) for b in arch_names[i+1:]]

summary = {
    'tool': 'EBI MAFFT REST (mafft)',
    'jobid': json.load(open('results/ebi_mafft_cpp_job.json'))['jobid'],
    'n_designed': len(des_names), 'n_archetypes': len(arch_names),
    'archetype_provenance': json.load(open('results/ebi_mafft_cpp_archetypes.json')),
    'max_identity_designed_vs_archetype': round(max(cross), 2),
    'mean_identity_designed_vs_archetype': round(sum(cross) / len(cross), 2),
    'mean_identity_within_designed': round(sum(within) / len(within), 2),
    'mean_identity_within_archetypes': round(sum(arch_within) / len(arch_within), 2),
    'n_designed_below_30pct_to_all_archetypes': sum(1 for d in des_names if per_designed[d]['nearest_identity'] < 30.0),
    'note': 'Ultra-short compositionally-biased peptides: alignment is descriptive, not phylogenetic.',
}
out = {'summary': summary, 'per_designed': per_designed}
json.dump(out, open('results/ebi_mafft_cpp_analysis.json', 'w'), indent=1)
print(json.dumps(summary, indent=1))
for d in des_names[:6]:
    p = per_designed[d]
    print(d, 'nearest:', p['nearest_archetype'], p['nearest_identity'])
