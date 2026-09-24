"""Structural verification of anchor-pocket preferences across HLA alleles.

For each structure: identify the peptide chain (shortest polymer, 8-12 aa),
the MHC alpha chain (longest chain excluding B2M), then compute which MHC
heavy-chain residues contact peptide position 2 (P2) and the C-terminal
position (P-omega) within 4.5 A (any heavy atom). Intersect contact sets
with the canonical B-pocket and F-pocket residue sets used in the 1DUZ
verification, and tabulate the amino-acid identity of every canonical
pocket residue in each allele. Output: results/pdb_pocket_verification_multi.json
"""
import json
from Bio.PDB import PDBParser

B_POCKET = [7, 9, 45, 63, 66, 67, 70, 99]
F_POCKET = [77, 80, 81, 84, 116, 123, 143, 146, 147]
CUTOFF = 4.5
THREE2ONE = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
    'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
    'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def analyze(pdb_id, path):
    st = PDBParser(QUIET=True).get_structure(pdb_id, path)[0]
    chains = {}
    for ch in st:
        res = [r for r in ch if r.id[0] == ' ' and r.resname in THREE2ONE]
        if res:
            chains[ch.id] = res
    # peptide = first chain of length 8..12 (multi-copy structures list
    # truncated copies too; the full-length peptide is the 8-12-mer)
    cands = [c for c, r in chains.items() if 8 <= len(r) <= 12]
    pep_id = min(cands, key=lambda c: -len(chains[c]))
    pep = chains[pep_id]
    pep_seq = ''.join(THREE2ONE[r.resname] for r in pep)

    # MHC alpha chain = the long chain (>150 aa) with most contacts to the
    # chosen peptide (picks the same copy in multi-copy asymmetric units)
    def n_contact(res_list, other):
        patoms = [a for r in other for a in r if a.element != 'H']
        n = 0
        for r in res_list:
            for a in r:
                if a.element == 'H':
                    continue
                if any((a - pa) <= CUTOFF for pa in patoms):
                    n += 1
                    break
        return n
    long_chains = [c for c, r in chains.items() if len(r) > 150 and c != pep_id]
    mhc_id = max(long_chains, key=lambda c: n_contact(chains[c], pep))
    mhc = {r.id[1]: r for r in chains[mhc_id]}
    mhc_seq_at = lambda i: THREE2ONE[mhc[i].resname] if i in mhc else None

    def contacts(pep_res):
        out = set()
        patoms = [a for a in pep_res if a.element != 'H']
        for rid, r in mhc.items():
            for a in r:
                if a.element == 'H':
                    continue
                for pa in patoms:
                    if (a - pa) <= CUTOFF:
                        out.add(rid)
                        break
                else:
                    continue
                break
        return sorted(out)

    p2 = contacts(pep[1])       # position 2
    pw = contacts(pep[-1])      # C-terminus
    return {
        'pdb': pdb_id,
        'chains': {'mhc': mhc_id, 'peptide': pep_id},
        'peptide_seq': pep_seq,
        'p2_contacts': p2,
        'pomega_contacts': pw,
        'p2_overlap_b': sorted(set(p2) & set(B_POCKET)),
        'pomega_overlap_f': sorted(set(pw) & set(F_POCKET)),
        'b_pocket_identities': {str(i): mhc_seq_at(i) for i in B_POCKET},
        'f_pocket_identities': {str(i): mhc_seq_at(i) for i in F_POCKET},
    }

out = {}
for pdb, path in [('1DUZ', 'data/raw/pdb/1duz.pdb'),
                  ('8RNI', 'data/raw/pdb/8rni.pdb'),
                  ('7OW3', 'data/raw/pdb/7ow3.pdb')]:
    out[pdb] = analyze(pdb, path)
json.dump(out, open('results/pdb_pocket_verification_multi.json', 'w'), indent=1)
for pdb, d in out.items():
    print(pdb, d['peptide_seq'], 'B-overlap', len(d['p2_overlap_b']),
          'F-overlap', len(d['pomega_overlap_f']))
print('B identities:')
for i in B_POCKET:
    print(i, [out[p]['b_pocket_identities'][str(i)] for p in out])
print('F identities:')
for i in F_POCKET:
    print(i, [out[p]['f_pocket_identities'][str(i)] for p in out])
