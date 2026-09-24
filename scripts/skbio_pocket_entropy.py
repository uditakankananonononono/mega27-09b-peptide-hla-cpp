"""scikit-bio Shannon entropy of pocket positions across the panel.

Question: is sequence variability across the 51-allele panel concentrated in
specific pocket positions, and where does residue 66 (the failed single-residue
epistasis-sign predictor, paper epistasis negative) sit in that ranking?
Uses skbio.diversity.alpha.shannon on per-position residue counts.
"""
import json
from collections import Counter
from skbio.diversity.alpha import shannon

B_POCKET = [7, 9, 24, 34, 45, 63, 66, 67, 70, 99]
F_POCKET = [74, 77, 80, 81, 84, 95, 97, 114, 116, 123, 143, 146, 147]
LEADER = 24

cov = [r["allele"] for r in json.load(open("results/per_allele_analysis.json"))["per_allele"]]
seqs, name, buf = {}, None, []
for line in open("data/raw/hla_prot.fasta"):
    line = line.strip()
    if line.startswith(">"):
        if name: seqs[name] = "".join(buf)
        p = line[1:].split(); name = p[1] if len(p) > 1 else None; buf = []
    elif line: buf.append(line)
if name: seqs[name] = "".join(buf)

panel = {}
for full, s in seqs.items():
    four = ":".join(full.split(":")[:2])
    if ("HLA-" + four) in cov and four not in panel:
        panel[four] = s

rows = []
for pocket, poss in (("B", B_POCKET), ("F", F_POCKET)):
    for p in poss:
        i = p - 1 + LEADER
        residues = [s[i] for s in panel.values() if len(s) > i]
        counts = list(Counter(residues).values())
        h = shannon(counts, base=2)
        rows.append({"pocket": pocket, "position": p, "n_alleles": len(residues),
                     "distinct_residues": len(counts), "shannon_bits": round(h, 3)})
rows.sort(key=lambda r: -r["shannon_bits"])
res66 = next(r for r in rows if r["position"] == 66)
out = {
    "tool": "scikit-bio (skbio.diversity.alpha.shannon, base 2)",
    "panel_alleles": len(panel),
    "n_positions": len(rows),
    "ranking_top5": rows[:5],
    "ranking_bottom5": rows[-5:],
    "res66_rank_of_23": [i + 1 for i, r in enumerate(rows) if r["position"] == 66][0],
    "res66": res66,
    "mean_entropy_B": round(sum(r["shannon_bits"] for r in rows if r["pocket"] == "B") / 10, 3),
    "mean_entropy_F": round(sum(r["shannon_bits"] for r in rows if r["pocket"] == "F") / 13, 3),
    "all": rows,
}
json.dump(out, open("results/hla_pocket_entropy_skbio.json", "w"), indent=1)
print("panel", len(panel), "| res66 rank", out["res66_rank_of_23"], "of 23, H =", res66["shannon_bits"],
      "| meanH B", out["mean_entropy_B"], "F", out["mean_entropy_F"])
print("top:", [(r["position"], r["shannon_bits"]) for r in rows[:3]])
