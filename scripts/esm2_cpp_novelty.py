"""ESM-2 protein language model embedding novelty analysis of the 18 CPPs.

Tool: facebookresearch ESM-2 (esm2_t6_8M_UR50D, 8M params, 320-d embeddings,
trained on UniRef50) via fair-esm. Each designed CPP, the 7 archetype CPPs
(from the MAFFT archetype fasta) and 300 natural 8-35-mer decoys
(UniProt reviewed short peptides) are embedded by mean-pooling the final
transformer layer. For each 9B-CPP we report the maximum cosine similarity
to any archetype CPP and its percentile rank against the decoy background -
a representation-space novelty score orthogonal to BLAST/phmmer alignment.
"""
import json, random, warnings
warnings.filterwarnings("ignore")
import torch, esm
import numpy as np

OUT = "results/cpp_esm2_novelty.json"
random.seed(7)

def load_cpp():
    cands = json.load(open("results/cpp_novel_candidates.json"))["named_candidates"]
    out = {}
    for c in cands:
        if isinstance(c, dict):
            out[c["name"]] = c["sequence"]
        else:
            out[c[0]] = c[1]
    return out

def load_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = ""
        elif name:
            seqs[name] += line.replace("-", "")
    return seqs

cpps = load_cpp()
arch = {k: v for k, v in load_fasta("results/ebi_mafft_cpp_alignment.fa").items()
        if k not in cpps}
decoy_pool = list(load_fasta("data/raw/uniprot_reviewed_len8_35.fasta").values())
decoys = random.sample(decoy_pool, min(300, len(decoy_pool)))

model, alphabet = esm.pretrained.esm2_t6_8M_UR50D()
model.eval()
bc = alphabet.get_batch_converter()

def embed(seqs, names):
    embs = {}
    B = 16
    for i in range(0, len(seqs), B):
        batch = [(names[j], seqs[j]) for j in range(i, min(i + B, len(seqs)))]
        _, _, toks = bc(batch)
        with torch.no_grad():
            res = model(toks, repr_layers=[6])
        rep = res["representations"][6]
        for k, (nm, sq) in enumerate(batch):
            embs[nm] = rep[k, 1:len(sq) + 1].mean(0).numpy()
    return embs

names = list(cpps) + list(arch) + [f"decoy{i}" for i in range(len(decoys))]
seqs = [cpps[n] for n in cpps] + [arch[n] for n in arch] + decoys
E = embed(seqs, names)

def cos(a, b):
    return float(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))

# background: each natural decoy's max cosine to the same 7 archetypes -
# identical neighbor count as the CPP query, so the percentile is fair
A = np.stack([E[a] for a in arch])
An = A / np.linalg.norm(A, axis=1, keepdims=True)
D = np.stack([E[f"decoy{i}"] for i in range(len(decoys))])
Dn = D / np.linalg.norm(D, axis=1, keepdims=True)
bg = np.sort((Dn @ An.T).max(axis=1))

rows = []
for n, sq in cpps.items():
    v = E[n]
    sims = {a: round(cos(v, E[a]), 4) for a in arch}
    best_arch = max(sims, key=sims.get)
    best = sims[best_arch]
    pct = float((bg < best).mean() * 100)
    rows.append({"name": n, "sequence": sq,
                 "nearest_archetype": best_arch, "max_cos_to_archetype": best,
                 "all_archetype_cos": sims,
                 "percentile_vs_decoy_maxcos": round(pct, 1)})

summary = {
    "model": "esm2_t6_8M_UR50D", "n_decoys": len(decoys),
    "archetypes": list(arch),
    "cpp_max_cos_to_archetype": {
        "mean": round(float(np.mean([r["max_cos_to_archetype"] for r in rows])), 4),
        "max": round(max(r["max_cos_to_archetype"] for r in rows), 4)},
    "decoy_background_maxcos": {"median": round(float(np.median(bg)), 4),
                                "p95": round(float(np.percentile(bg, 95)), 4)},
    "n_cpp_below_decoy_p95": sum(1 for r in rows if r["max_cos_to_archetype"] < np.percentile(bg, 95)),
}
json.dump({"rows": rows, "summary": summary}, open(OUT, "w"), indent=1)
print(json.dumps(summary, indent=1))
for r in sorted(rows, key=lambda r: r["max_cos_to_archetype"]):
    print(r["name"], "->", r["nearest_archetype"], r["max_cos_to_archetype"],
          "pct", r["percentile_vs_decoy_maxcos"])
