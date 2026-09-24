"""EMBOSS antigenic (EBI REST) - Kolaskar-Tongaonkar per-residue physicochemical profiles of the 18 CPPs.

EMBOSS antigenic predicts antigenic determinants (Kolaskar & Tongaonkar 1990,
semi-empirical residue antigenicity scale). Used here as an orthogonal
immunogenicity screen on the designed CPPs - a peptide that is both
cell-penetrating and strongly antigenic is a worse delivery vehicle.

Usage:
  python3 scripts/emboss_pepinfo_screen.py submit
  python3 scripts/emboss_pepinfo_screen.py poll
"""
import json, os, sys, urllib.request, urllib.parse

BASE = "https://www.ebi.ac.uk/Tools/services/rest/emboss_pepinfo"
EMAIL = "uditakankana@gmail.com"
CACHE = "/tmp/pepinfo_raw"
os.makedirs(CACHE, exist_ok=True)
JOBS_PATH = "results/cpp_emboss_pepinfo_jobs.json"
OUT_PATH = "results/cpp_emboss_pepinfo.json"

def get(url, timeout=30):
    req = urllib.request.Request(url, headers={"User-Agent": "mega27-09b"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

def post(url, data, timeout=30):
    body = urllib.parse.urlencode(data).encode()
    req = urllib.request.Request(url, data=body, headers={"User-Agent": "mega27-09b"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace").strip()

def load_seqs():
    cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
    return {c["name"]: c["sequence"] if isinstance(c, dict) else (c[0], c[1])
            for c in cands} if isinstance(cands[0], dict) else dict(cands)

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
    print("submitted:", len(jobs))

def poll():
    jobs = json.load(open(JOBS_PATH))
    results = json.load(open(OUT_PATH)) if os.path.exists(OUT_PATH) else {}
    done = queued = err = 0
    for name, j in sorted(jobs.items()):
        if name in results:
            done += 1
            continue
        jid = j["jobid"]
        try:
            st = get(f"{BASE}/status/{jid}").strip()
        except Exception as e:
            print(name, "poll error", e)
            err += 1
            continue
        if st == "FINISHED":
            try:
                out = get(f"{BASE}/result/{jid}/out")
                open(f"{CACHE}/{name}.out", "w").write(out)
                results[name] = {"jobid": jid, "raw_file": f"{CACHE}/{name}.out"}
                done += 1
                print(name, "DONE", len(out), "bytes", flush=True)
            except Exception as e:
                print(name, "fetch error", e)
                err += 1
        else:
            queued += 1
            print(name, st, flush=True)
    json.dump(results, open(OUT_PATH, "w"), indent=1)
    print(f"done={done} queued={queued} err={err} -> {OUT_PATH}")

if __name__ == "__main__":
    {"submit": submit, "poll": poll}[sys.argv[1]]()
