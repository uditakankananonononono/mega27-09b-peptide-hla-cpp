"""EMBOSS pepstats (EBI REST) physicochemical screen of the 18 CPP candidates.

Tool 40 in the external-tools audit: EBI EMBOSS pepstats computes molecular
weight, isoelectric point, charge, extinction coefficients and amino-acid
composition per candidate. REST API:
https://www.ebi.ac.uk/Tools/services/rest/emboss_pepstats

Usage:
  python3 scripts/emboss_pepstats_screen.py submit   # submit all jobs
  python3 scripts/emboss_pepstats_screen.py poll     # poll + collect finished
"""
import json, os, sys, time, urllib.request, urllib.parse

BASE = "https://www.ebi.ac.uk/Tools/services/rest/emboss_pepstats"
EMAIL = "uditakankana@gmail.com"
CACHE = "/tmp/pepstats_raw"
os.makedirs(CACHE, exist_ok=True)
JOBS_PATH = "results/cpp_emboss_pepstats_jobs.json"
OUT_PATH = "results/cpp_emboss_pepstats.json"

def get(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "mega27-09b"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

def post(url, data, timeout=30):
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, data=body,
                                 headers={"User-Agent": "mega27-09b"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace").strip()

def load_seqs():
    cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
    seqs = {}
    for c in cands:
        if isinstance(c, dict):
            seqs[c["name"]] = c["sequence"]
        else:
            seqs[c[0]] = c[1]
    return seqs

def submit():
    seqs = load_seqs()
    jobs = json.load(open(JOBS_PATH)) if os.path.exists(JOBS_PATH) else {}
    for name, seq in sorted(seqs.items()):
        if name in jobs:
            continue
        jobs[name] = {"jobid": post(BASE + "/run", {"email": EMAIL, "sequence": seq}),
                      "sequence": seq}
        print(name, "->", jobs[name]["jobid"], flush=True)
        json.dump(jobs, open(JOBS_PATH, "w"), indent=1)
    json.dump(jobs, open(JOBS_PATH, "w"), indent=1)
    print("submitted:", len(jobs))

def poll():
    jobs = json.load(open(JOBS_PATH))
    results = json.load(open(OUT_PATH)) if os.path.exists(OUT_PATH) else {}
    done = queued = err = 0
    for name, j in sorted(jobs.items()):
        cf = os.path.join(CACHE, name + ".out")
        if os.path.exists(cf):
            raw = open(cf).read()
        else:
            jid = j["jobid"] if isinstance(j, dict) else j
            st = get(BASE + "/status/" + jid).strip()
            if st != "FINISHED":
                print(name, st, flush=True)
                if st in ("ERROR", "FAILURE", "NOT_FOUND"):
                    err += 1
                else:
                    queued += 1
                continue
            raw = get(BASE + "/result/" + jid + "/out")
            open(cf, "w").write(raw)
        results[name] = {"sequence": j["sequence"] if isinstance(j, dict) else load_seqs()[name],
                         "raw": raw}
        done += 1
        print(name, "DONE", len(raw), "bytes", flush=True)
    json.dump(results, open(OUT_PATH, "w"), indent=1)
    print(f"done={done} queued={queued} err={err} -> {OUT_PATH}")

if __name__ == "__main__":
    (submit if (len(sys.argv) > 1 and sys.argv[1] == "submit") else poll)()
