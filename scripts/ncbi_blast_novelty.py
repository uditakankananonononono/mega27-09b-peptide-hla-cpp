"""NCBI BLAST protein-level novelty screen of the 18 designed CPPs (tool: NCBI BLAST URLAPI).

Polls the 18 queued blastp (Swiss-Prot) jobs in results/cpp_ncbi_blast_jobs.json,
fetches tabular alignments for jobs that are READY, and writes per-peptide top
hits plus a novelty summary to results/cpp_ncbi_blast.json. Swiss-Prot subject
accessions of genuine hits are accession-level datasets under the program rule.

Run: python3 scripts/ncbi_blast_novelty.py   (requires network access)
"""
import json, time, urllib.parse, urllib.request

BASE = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi'


def status(rid):
    url = f'{BASE}?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}'
    txt = urllib.request.urlopen(url, timeout=20).read().decode()
    if 'Status=READY' in txt:
        return 'READY'
    if 'Status=WAITING' in txt:
        return 'WAITING'
    return 'UNKNOWN'


def fetch_tabular(rid):
    url = f'{BASE}?CMD=Get&FORMAT_TYPE=Text&ALIGNMENT_VIEW=Tabular&RID={rid}'
    return urllib.request.urlopen(url, timeout=60).read().decode()


def parse_tabular(text):
    hits = []
    for line in text.splitlines():
        if line.startswith('#') or not line.strip():
            continue
        f = line.split('\t')
        if len(f) < 12:
            continue
        hits.append({'subject_id': f[1], 'pct_identity': float(f[2]),
                     'align_len': int(f[3]), 'evalue': float(f[10]),
                     'bit_score': float(f[11])})
    return hits


def main():
    jobs = json.load(open('results/cpp_ncbi_blast_jobs.json'))
    try:
        out = json.load(open('results/cpp_ncbi_blast.json'))
    except FileNotFoundError:
        out = {'tool': 'NCBI BLAST URLAPI (blastp, Swiss-Prot, short-query)',
               'tool_url': 'https://blast.ncbi.nlm.nih.gov/Blast.cgi',
               'jobs': jobs, 'results': {}}
    for name, rid in sorted(jobs.items()):
        if name in out['results'] and out['results'][name].get('status') == 'READY':
            continue
        st = status(rid)
        if st != 'READY':
            out['results'][name] = {'status': st, 'rid': rid}
            print(name, st)
            continue
        hits = parse_tabular(fetch_tabular(rid))
        out['results'][name] = {'status': 'READY', 'rid': rid, 'n_hits': len(hits),
                                'top_hits': hits[:5]}
        print(name, 'READY', len(hits), 'hits')
        time.sleep(1)  # NCBI usage policy: <=1 request per 2s is safer; keep polite
    n_ready = sum(1 for r in out['results'].values() if r.get('status') == 'READY')
    out['summary'] = {'n_ready': n_ready, 'n_total': len(jobs)}
    json.dump(out, open('results/cpp_ncbi_blast.json', 'w'), indent=1)


if __name__ == '__main__':
    main()
