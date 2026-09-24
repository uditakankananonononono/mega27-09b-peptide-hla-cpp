"""PDB-REDO QC of the full per-allele structure set: fetch data.json per PDB
entry, record resolution and R-free, flag any that fail the paper's QC bar
(resolution <= 3.5 A, R-free <= 0.35). Each PDB-REDO entry is an accession-level
dataset under the program counting rule.
"""
import json, os, time, urllib.request

rows = json.load(open("results/pdb_per_allele_structures.json"))
pdbs = sorted({r["pdb"] for r in rows})
out = {}
for pdb in pdbs:
    fn = f"data/raw/pdb_structures/{pdb}.redo.json"
    if not os.path.exists(fn):
        url = f"https://pdb-redo.eu/db/{pdb}/data.json"
        for i in range(3):
            try:
                urllib.request.urlretrieve(url, fn)
                break
            except Exception as e:
                print("retry", pdb, i, str(e)[:80], flush=True)
                time.sleep(2 * (i + 1))
        time.sleep(0.3)
    try:
        d = json.load(open(fn))
        props = d.get("properties", {})
        out[pdb] = {
            "resolution": props.get("RESOLUTION"),
            "r_free": props.get("RFFIN") or props.get("RFREE"),
            "r_work": props.get("RFIN") or props.get("R"),
        }
        print(pdb, out[pdb]["resolution"], out[pdb]["r_free"], flush=True)
    except Exception as e:
        out[pdb] = {"error": str(e)[:120]}
        print("ERR", pdb, str(e)[:100], flush=True)

ok = [v for v in out.values() if "error" not in v and v.get("resolution")]
fails = {p: v for p, v in out.items() if "error" not in v and v.get("resolution")
         and (float(v["resolution"]) > 3.5 or (v.get("r_free") and float(v["r_free"]) > 0.35))}
summary = {"n_qc": len(out), "n_ok": len(ok), "n_qc_fail": len(fails), "fails": fails}
json.dump({"per_structure": out, "summary": summary}, open("results/structure_qc_all.json", "w"), indent=1)
print("SUMMARY", summary["n_qc"], summary["n_ok"], "fails:", list(fails))
