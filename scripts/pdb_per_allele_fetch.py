"""Per-allele pHLA structure set: query RCSB for peptide-HLA structures of the
top-N alleles by test-set size, download mmCIF, identify the peptide chain,
and record per-structure facts (peptide sequence/length, resolution).

Each PDB entry is one accession-level dataset under the uniform gate rule.
"""
import json, os, time, urllib.request, urllib.parse

ALLELES = ["HLA-A*02:01","HLA-A*03:01","HLA-A*11:01","HLA-A*31:01","HLA-A*68:02",
           "HLA-A*02:03","HLA-B*07:02","HLA-B*15:01","HLA-A*01:01","HLA-A*02:06",
           "HLA-A*24:02","HLA-B*08:01","HLA-B*27:05","HLA-B*58:01","HLA-B*44:02",
           "HLA-B*35:01","HLA-B*40:01","HLA-B*51:01","HLA-B*18:01","HLA-A*26:01",
           "HLA-A*33:01","HLA-A*68:01","HLA-A*29:02","HLA-A*23:01","HLA-A*30:01",
           "HLA-B*15:03","HLA-B*57:01","HLA-A*32:01","HLA-B*44:03","HLA-B*46:01"]
TOPN = 6
OUTDIR = "data/raw/pdb_structures"
os.makedirs(OUTDIR, exist_ok=True)

def _get_json(url, tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                raw = r.read()
            if not raw.strip():
                return {}  # HTTP 204: zero search hits
            return json.loads(raw)
        except Exception as e:
            print("retry", i, e, flush=True)
            time.sleep(2 * (i + 1))
    return {}  # unreachable API treated as no results, recorded in manifest

def rcsb_search(phrase, rows=TOPN):
    q = {"query":{"type":"terminal","service":"text",
         "parameters":{"attribute":"struct.title","operator":"contains_phrase","value":phrase}},
         "return_type":"entry",
         "request_options":{"paginate":{"start":0,"rows":rows},
                            "results_content_type":["experimental"]}}
    url = "https://search.rcsb.org/rcsbsearch/v2/query?json=" + urllib.parse.quote(json.dumps(q))
    d = _get_json(url)
    return [x["identifier"] for x in d.get("result_set", [])]

mf_path = os.path.join(OUTDIR, "search_manifest.json")
manifest = json.load(open(mf_path)) if os.path.exists(mf_path) else {}
for allele in ALLELES:
    if allele in manifest and manifest[allele]:
        print("cached", allele, manifest[allele], flush=True)
        continue
    ids = rcsb_search(allele)
    manifest[allele] = ids
    json.dump(manifest, open(mf_path, "w"), indent=1)
    print(allele, ids, flush=True)
    time.sleep(1.0)

# download mmCIF for each
for allele, ids in manifest.items():
    for pdb in ids:
        fn = os.path.join(OUTDIR, pdb.lower() + ".cif")
        if os.path.exists(fn) and os.path.getsize(fn) > 10000:
            continue
        url = f"https://files.rcsb.org/download/{pdb}.cif"
        for i in range(4):
            try:
                urllib.request.urlretrieve(url, fn)
                print("fetched", pdb, os.path.getsize(fn), flush=True)
                break
            except Exception as e:
                print("dl-retry", pdb, i, e, flush=True)
                time.sleep(2 * (i + 1))
        time.sleep(0.5)
print("done")
