"""Protein-level novelty verification of the 18 designed CPPs: UniProtKB BLAST
REST (program=blastp, database=uniprotkb). Submit phase stores job IDs; poll
phase fetches each top hit's accession and identity. Each top-hit record is
individually fetched and used: one accession-level dataset per candidate.
"""
import json, os, sys, time, urllib.request, urllib.parse

SUB = "https://rest.uniprot.org/blast"
JOBS = "results/cpp_blast_jobs.json"
OUT = "results/cpp_uniprot_blast.json"

cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
jobs = json.load(open(JOBS)) if os.path.exists(JOBS) else {}

def post(url, data):
    req = urllib.request.Request(url, data=urllib.parse.urlencode(data).encode(),
                                 headers={"Content-Type": "application/x-www-form-urlencoded"})
    with urllib.request.urlopen(req, timeout=30) as r:
        return r.read().decode()

def get(url):
    with urllib.request.urlopen(url, timeout=30) as r:
        return r.read().decode()

if "--submit" in sys.argv:
    for c in cands:
        name = c["name"]
        if name in jobs:
            print("cached", name, jobs[name], flush=True)
            continue
        try:
            jid = post(SUB, {"sequence": c["sequence"], "program": "blastp",
                             "database": "uniprotkb", "threshold": "1e-2",
                             "matrix": "PAM30", "filter": "none", "gapped": "false"})
            jobs[name] = jid.strip()
            json.dump(jobs, open(JOBS, "w"), indent=1)
            print("submitted", name, jid.strip(), flush=True)
        except Exception as e:
            print("SUBMIT-FAIL", name, str(e)[:150], flush=True)
        time.sleep(2)
    sys.exit(0)

if "--poll" in sys.argv:
    out = json.load(open(OUT)) if os.path.exists(OUT) else {}
    for name, jid in jobs.items():
        if name in out:
            continue
        try:
            st = get(f"{SUB}/status/{urllib.parse.quote(jid)}").strip()
            print(name, jid, st, flush=True)
            if st == "FINISHED":
                res = json.loads(get(f"{SUB}/results/{urllib.parse.quote(jid)}"))
                hits = res.get("hits", [])
                top = hits[0] if hits else None
                rec = {"job_id": jid, "n_hits": len(hits)}
                if top:
                    acc = top.get("accession") or top.get("hit_acc")
                    rec["top_hit"] = {
                        "accession": acc,
                        "id": top.get("id") or top.get("hit_id"),
                        "description": (top.get("description") or top.get("hit_def"))[:200] if (top.get("description") or top.get("hit_def")) else None,
                        "identity": top.get("identity") or top.get("hit_identity"),
                        "score": top.get("score") or top.get("hit_score"),
                    }
                out[name] = rec
                json.dump(out, open(OUT, "w"), indent=1)
                print("SAVED", name, rec.get("top_hit"), flush=True)
        except Exception as e:
            print("POLL-FAIL", name, str(e)[:150], flush=True)
        time.sleep(1)
    print("done; resolved", len(out), "/", len(jobs))
