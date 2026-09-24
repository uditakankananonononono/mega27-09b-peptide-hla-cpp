"""Analyze the per-allele pHLA structure set: for each mmCIF, identify the
peptide chain (8-12 residues, not the ~275-aa heavy chain or ~99-aa B2M),
extract sequence, resolution, and per-residue relative solvent accessibility
(freesasa). Tests whether P2 and P-Omega anchor burial generalizes beyond the
original 3-structure set.
"""
import json, os, glob, warnings, statistics as stats
warnings.filterwarnings("ignore")
from Bio.PDB import MMCIFParser, MMCIF2Dict, PDBIO
import freesasa, numpy as np

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E',
       'GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F',
       'PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
ALLELE_OF = {
 '7t5m':'HLA-A*02:01','9nmu':'HLA-A*02:01','8fu4':'HLA-A*02:01','9nmv':'HLA-A*02:01','9nmy':'HLA-A*02:01','9ytd':'HLA-A*02:01',
 '8rni':'HLA-A*03:01','8vjz':'HLA-A*03:01','8dvg':'HLA-A*03:01','7l1c':'HLA-A*03:01','7stf':'HLA-A*03:01','8vcl':'HLA-A*03:01',
 '5wjn':'HLA-A*11:01','5wjl':'HLA-A*11:01','5wkf':'HLA-A*11:01','5wkh':'HLA-A*11:01','9wpd':'HLA-A*11:01','7ow3':'HLA-A*11:01',
 '3ox8':'HLA-A*02:03','6uj7':'HLA-B*07:02','6uj8':'HLA-B*07:02','7kgu':'HLA-B*07:02','7s7e':'HLA-B*07:02','7s7f':'HLA-B*07:02','6avf':'HLA-B*07:02',
 '6uzp':'HLA-B*15:01','6uzq':'HLA-B*15:01','6uzs':'HLA-B*15:01','6vb3':'HLA-B*15:01','8elg':'HLA-B*15:01','8elh':'HLA-B*15:01',
 '6mpp':'HLA-A*01:01','9yir':'HLA-A*01:01','9z50':'HLA-A*01:01','3oxr':'HLA-A*02:06'}
MAXASA = {'A':121,'R':265,'N':187,'D':187,'C':148,'Q':214,'E':214,'G':97,'H':216,
          'I':195,'L':191,'K':230,'M':203,'F':228,'P':154,'S':143,'T':163,'W':264,'Y':255,'V':165}

def chain_seq(chain):
    seq = []
    for res in chain:
        rn = res.get_resname().strip()
        if res.id[0] != ' ':
            continue
        if rn in AA3:
            seq.append(AA3[rn])
        elif rn not in ('HOH',):
            return None
    return ''.join(seq) if seq else None

parser = MMCIFParser(QUIET=True)
rows = []
for fn in sorted(glob.glob('data/raw/pdb_structures/*.cif')):
    pdb = os.path.basename(fn)[:4]
    allele = ALLELE_OF.get(pdb)
    try:
        d = MMCIF2Dict.MMCIF2Dict(fn)
        reso = d.get('_refine.ls_d_res_high', [None])[0]
        reso = float(reso) if reso not in (None, '.', '?') else None
        st = parser.get_structure(pdb, fn)
        model = st[0]
        chains = []
        for ch in model:
            s = chain_seq(ch)
            if s:
                chains.append((ch.id, len(s), s))
        pep = [c for c in chains if 8 <= c[1] <= 12]
        if not pep:
            rows.append({'pdb': pdb, 'allele': allele, 'resolution': reso,
                         'peptide': None, 'note': 'no 8-12mer peptide chain',
                         'chains': [(c[0], c[1]) for c in chains]})
            print(pdb, allele, reso, 'NO PEPTIDE', [(c[0], c[1]) for c in chains][:6], flush=True)
            continue
        cid, plen, pseq = min(pep, key=lambda c: c[1])
        tmp = f'/tmp/{pdb}.pdb'
        io = PDBIO(); io.set_structure(model); io.save(tmp)
        r = freesasa.calc(freesasa.Structure(tmp))
        ra = r.residueAreas()
        vals = list(ra.get(cid, {}).values())
        total = [v.total for v in vals]
        rsa = [round(t / MAXASA.get(pseq[i], 200), 3) for i, t in enumerate(total[:plen])]
        rows.append({'pdb': pdb, 'allele': allele, 'resolution': reso,
                     'peptide': pseq, 'pep_chain': cid, 'pep_len': plen, 'rsa': rsa,
                     'p2_rsa': rsa[1] if plen >= 2 else None,
                     'pomega_rsa': rsa[-1] if rsa else None,
                     'min_rsa_pos': int(np.argmin(rsa)) + 1 if rsa else None,
                     'chains': [(c[0], c[1]) for c in chains]})
        print(pdb, allele, reso, pseq, 'P2', rsa[1] if plen >= 2 else None, 'PO', rsa[-1], flush=True)
    except Exception as e:
        rows.append({'pdb': pdb, 'allele': allele, 'error': str(e)[:200]})
        print('ERR', pdb, str(e)[:120], flush=True)

json.dump(rows, open('results/pdb_per_allele_structures.json', 'w'), indent=1)
ok = [r for r in rows if r.get('rsa')]
if not ok:
    print('NO STRUCTURES ANALYZED'); raise SystemExit(1)
n_ok = len(ok)
p2_top3 = sum(1 for r in ok if 1 in sorted(range(len(r['rsa'])), key=lambda i: r['rsa'][i])[:3])
po_top3 = sum(1 for r in ok if (len(r['rsa'])-1) in sorted(range(len(r['rsa'])), key=lambda i: r['rsa'][i])[:3])
internal = [x for r in ok for i, x in enumerate(r['rsa']) if 2 <= i < len(r['rsa'])-1]
summary = {'n_structures': len(rows), 'n_with_peptide': n_ok,
           'n_no_peptide': sum(1 for r in rows if r.get('peptide') is None and 'error' not in r),
           'n_errors': sum(1 for r in rows if 'error' in r),
           'p2_top3_buried': p2_top3, 'pomega_top3_buried': po_top3,
           'mean_p2_rsa': round(stats.mean(r['p2_rsa'] for r in ok), 4),
           'mean_pomega_rsa': round(stats.mean(r['pomega_rsa'] for r in ok), 4),
           'mean_internal_rsa': round(stats.mean(internal), 4),
           'alleles': sorted({r['allele'] for r in ok})}
json.dump(summary, open('results/pdb_per_allele_summary.json', 'w'), indent=1)
print(json.dumps(summary, indent=1))
