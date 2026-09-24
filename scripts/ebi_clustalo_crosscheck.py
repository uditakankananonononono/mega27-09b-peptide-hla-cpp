"""Submit/poll EBI Clustal Omega (REST) for the 18 designed CPPs + 7 CPPsite archetypes.

Cross-engine check on the MAFFT family-level relatedness result: does an
independent MSA engine (Clustal Omega) reproduce the low identity of the
designed CPPs to known CPP archetypes? Usage: submit | poll
"""
import json, os, sys, time, urllib.parse, urllib.request

BASE = 'https://www.ebi.ac.uk/Tools/services/rest/clustalo'
EMAIL = 'uditakankana@gmail.com'
FASTA = 'results/ebi_mafft_cpp_alignment.fa'
JOB = 'results/cpp_ebi_clustalo_job.json'
RAW = '/tmp/clustalo_raw'
os.makedirs(RAW, exist_ok=True)

def submit():
    seq = open(FASTA).read()
    data = urllib.parse.urlencode({'email': EMAIL, 'stype': 'protein',
                                   'sequence': seq}).encode()
    req = urllib.request.Request(BASE + '/run', data=data)
    jobid = urllib.request.urlopen(req, timeout=60).read().decode().strip()
    json.dump({'jobid': jobid, 'submitted': time.strftime('%Y-%m-%dT%H:%M:%S')},
              open(JOB, 'w'), indent=1)
    print('submitted', jobid)

def poll():
    j = json.load(open(JOB)); jobid = j['jobid']
    st = urllib.request.urlopen(f'{BASE}/status/{jobid}', timeout=30).read().decode().strip()
    print(jobid, st)
    if st == 'FINISHED' and not j.get('done'):
        for rt in ('aln-clustal_num', 'pim', 'phylotree'):
            try:
                r = urllib.request.urlopen(f'{BASE}/result/{jobid}/{rt}', timeout=60).read()
                open(f'{RAW}/{rt}.out', 'wb').write(r)
                print(rt, len(r), 'bytes')
            except Exception as e:
                print(rt, 'ERR', e)
        j['done'] = True
        json.dump(j, open(JOB, 'w'), indent=1)

if __name__ == '__main__':
    {'submit': submit, 'poll': poll}[sys.argv[1]]()
