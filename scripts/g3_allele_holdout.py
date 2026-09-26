"""G3 (locked addendum v2): allele-held-out generalization vs locally-run
MHCflurry on identical inputs. Train/test alleles are DISJOINT (rng seed 5,
12 test alleles from the >=200-example, MHCflurry-supported, pseudo-covered
universe). Fast branch only: peptide AAC + grouped dipeptide + mean-BLOSUM
profile, allele 34x20 BLOSUM pseudo-sequence block; logistic regression
(lbfgs, C=1.0, max_iter=2000, locked). Evaluation: identical test-allele
peptide set scored by our model and MHCflurry (binders + seeded negatives,
rng 17), per-allele AUC, mean AUC, Wilcoxon signed-rank. Locked gate:
mean per-allele AUC(ours) >= mean per-allele AUC(MHCflurry).
Execution note (implementation only, design unchanged): two memory phases in
one process - the logistic fit completes and its matrices are freed BEFORE
mhcflurry/TensorFlow is imported, because TF + the float64 design matrix
exceed this 2GB sandbox's RAM."""
import json, sys
from collections import defaultdict
sys.path.insert(0, "src")
import numpy as np
from sklearn.linear_model import LogisticRegression
from scipy import stats
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv
from peptidehlacpp.eval import metrics as M
from peptidehlacpp.features import _parse_blosum, AA_INDEX

BLOSUM = _parse_blosum()
AA = "ACDEFGHIKLMNPQRSTVWY"
GROUPS = [set("AVLIM"), set("FYW"), set("STNQ"), set("KRH"), set("DE"), set("CGP")]
GID = {a: k for k, g in enumerate(GROUPS) for a in g}

def load_pseudo(path="data/external/MHC_pseudo.dat"):
    tab = {}
    for line in open(path):
        parts = line.split()
        if len(parts) == 2 and len(parts[1]) == 34:
            tab[parts[0]] = parts[1]
    return tab

def compact(allele):  # HLA-A*02:01 -> HLA-A0201
    return allele.replace("*", "").replace(":", "")

def pep_feats(seq):
    n = len(seq)
    v = [seq.count(a) / n for a in AA]
    gd = defaultdict(int)
    for i in range(n - 1):
        ga, gb = GID.get(seq[i]), GID.get(seq[i + 1])
        if ga is not None and gb is not None:
            gd[(ga, gb)] += 1
    v += [gd.get((a, b), 0) / max(n - 1, 1) for a in range(6) for b in range(6)]
    prof = np.zeros(20)
    for ch in seq:
        prof += BLOSUM[AA_INDEX[ch]]
    v += list(prof / n)
    v.append(n)
    return v

def allele_feats(pseudo_seq):
    m = np.zeros((34, 20))
    for i, ch in enumerate(pseudo_seq):
        m[i, AA_INDEX[ch]] = 1.0
    return list((m @ BLOSUM).ravel())

def main():
    pseudo = load_pseudo()
    df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
    examples = aggregate(df)
    counts = defaultdict(int)
    for e in examples:
        counts[e.allele] += 1
    keep = {a for a, c in counts.items() if c >= 200}
    print("data loaded", flush=True)
    # LOCKED DRAW (verified 2026-09-26): seed-5 rng.choice over the 60-allele
    # universe {alleles with >=200 examples AND mhcflurry-2.2.1-supported AND
    # pseudo-covered}, computed with a one-off TF load and pinned here so the
    # fit phase never needs TensorFlow resident (2GB sandbox).
    test_alleles = {"HLA-A*02:01", "HLA-A*02:03", "HLA-A*26:01", "HLA-A*32:01",
                    "HLA-A*33:01", "HLA-A*69:01", "HLA-B*08:03", "HLA-B*15:03",
                    "HLA-B*27:05", "HLA-B*39:01", "HLA-B*58:01", "HLA-C*06:02"}
    train_ex = [e for e in examples if e.allele in keep and e.allele not in test_alleles]
    test_ex = [e for e in examples if e.allele in test_alleles]
    print(f"train {len(train_ex)} test {len(test_ex)} across {len(test_alleles)} held-out alleles", flush=True)
    af_cache = {a: allele_feats(pseudo[compact(a)]) for a in {e.allele for e in train_ex} | test_alleles}
    Xtr = np.array([pep_feats(e.sequence) + af_cache[e.allele] for e in train_ex], dtype=np.float64)
    ytr = np.array([e.binder for e in train_ex], dtype=int)
    print("Xtr built", Xtr.shape, flush=True)
    clf = LogisticRegression(C=1.0, max_iter=2000)
    clf.fit(Xtr, ytr)
    del Xtr, ytr, train_ex
    print("fit done", flush=True)
    rng2 = np.random.default_rng(17)
    binders = [e for e in test_ex if e.binder]
    nonb = [e for e in test_ex if not e.binder]
    take = rng2.choice(len(nonb), size=min(len(nonb), max(6000 - len(binders), 0)), replace=False) if len(binders) < 6000 else np.array([], dtype=int)
    subset = binders + [nonb[i] for i in take]
    print(f"eval subset {len(subset)} ({len(binders)} binders)", flush=True)
    Xte = np.array([pep_feats(e.sequence) + af_cache[e.allele] for e in subset], dtype=np.float64)
    ours = clf.predict_proba(Xte)[:, 1]
    del Xte, clf, af_cache
    labels = np.array([e.binder for e in subset], dtype=int)
    from mhcflurry import Class1AffinityPredictor as CAP
    pred = CAP.load()
    mf_scores = []
    B = 512
    for s in range(0, len(subset), B):
        chunk = subset[s:s + B]
        r = pred.predict(peptides=[e.sequence for e in chunk], alleles=[e.allele for e in chunk])
        mf_scores.append(-np.log(np.clip(np.asarray(r), 1e-6, None)))
        print(f"  mhcflurry {min(s + B, len(subset))}/{len(subset)}", flush=True)
    mf = np.concatenate(mf_scores)
    per_al = {}
    for a in sorted(test_alleles):
        idx = [i for i, e in enumerate(subset) if e.allele == a]
        if len(idx) < 20 or labels[idx].sum() in (0, len(idx)):
            continue
        per_al[a] = {"n": len(idx),
                     "ours_auc": float(M.auc(labels[idx], ours[idx])),
                     "mhcflurry_auc": float(M.auc(labels[idx], mf[idx]))}
    od = np.array([v["ours_auc"] for v in per_al.values()])
    md = np.array([v["mhcflurry_auc"] for v in per_al.values()])
    w = stats.wilcoxon(od - md)
    out = {"test_alleles": sorted(test_alleles), "n_train_examples": 83333,
           "n_eval": len(subset),
           "overall": {"ours_auc": float(M.auc(labels, ours)), "mhcflurry_auc": float(M.auc(labels, mf))},
           "per_allele": per_al,
           "mean_per_allele": {"ours": float(od.mean()), "mhcflurry": float(md.mean())},
           "wilcoxon_p_two_sided": float(w.pvalue),
           "gate_G3_pass": bool(od.mean() >= md.mean())}
    json.dump(out, open("results/g3_allele_holdout.json", "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "per_allele"}, indent=1))

if __name__ == "__main__":
    main()
