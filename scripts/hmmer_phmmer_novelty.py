"""EBI HMMER phmmer (REST v1) novelty screen of the 18 designed CPPs.

Orthogonal novelty engine to NCBI BLAST: profile-HMM based sequence search
of each designed CPP against Swiss-Prot. submit | poll, resumable.
"""
import json, os, sys, time, urllib.request

BASE = 'https://www.ebi.ac.uk/Tools/hmmer/api/v1/search'
JOBS = 'results/cpp_hmmer_phmmer_jobs.json'
OUT = 'results/cpp_hmmer_phmmer.json'
RAW = '/tmp/hmmer_raw'
os.makedirs(RAW, exist_ok=True)

def cands():
    out = []
    for c in json.load(open('results/cpp_novel_candidates.json'))['named_candidates']:
        out.append((c['name'], c['sequence']) if isinstance(c, dict) else (c[0], c[1]))
    return out

def submit():
    jobs = json.load(open(JOBS)) if os.path.exists(JOBS) else {}
    for name, seq in cands():
        if name in jobs: continue
        body = json.dumps({'input': f'>{name}\n{seq}', 'database': 'swissprot'}).encode()
        req = urllib.request.Request(f'{BASE}/phmmer', data=body,
                                     headers={'Content-Type': 'application/json',
                                              'Accept': 'application/json'})
        d = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
        jobs[name] = d['id']
        print('submitted', name, d['id'])
        time.sleep(1)
    json.dump(jobs, open(JOBS, 'w'), indent=1)

def poll():
    jobs = json.load(open(JOBS))
    res = json.load(open(OUT)) if os.path.exists(OUT) else {}
    done = 0
    for name, jid in jobs.items():
        if name in res: done += 1; continue
        req = urllib.request.Request(f'{BASE}/{jid}', headers={'Accept': 'application/json'})
        try:
            d = json.loads(urllib.request.urlopen(req, timeout=60).read().decode())
        except Exception as e:
            print(name, 'ERR', str(e)[:80]); continue
        if 'results' not in d:
            print(name, 'RUNNING'); continue
        open(f'{RAW}/{name}.json', 'w').write(json.dumps(d))
        hits = d['results'].get('hits', [])
        res[name] = {'n_hits': len(hits),
                     'top': [{'acc': h.get('acc'), 'name': h.get('name'),
                              'evalue': h.get('evalue'),
                              'desc': (h.get('desc') or '')[:80]} for h in hits[:5]]}
        done += 1
        print(name, 'DONE hits', len(hits))
    json.dump(res, open(OUT, 'w'), indent=1)
    print(f'done={done}/{len(jobs)}')

if __name__ == '__main__':
    {'submit': submit, 'poll': poll}[sys.argv[1]]()
