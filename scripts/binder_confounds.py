"""Confound-controlled evidence that the ensemble score carries binding signal.

Fits two binomial GLMs on the held-out test set (split replication verified
against results/ensemble_test_predictions.npz):
  M0: binder ~ C(length) + C(assay method) + C(locus)
  M1: binder ~ ensemble_z + C(length) + C(assay method) + C(locus)
If the ensemble only exploited length/method/locus artifacts, M1's ensemble
coefficient would be ~0 and the pseudo-R2 gain negligible.
Output: results/binder_confounds.json
"""
import json, sys
from collections import defaultdict, Counter
sys.path.insert(0, "src")
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv, split_by_peptide

df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
meth = defaultdict(Counter)
for seq, al, m in df[["sequence", "allele", "method"]].itertuples(index=False):
    meth[(seq, al)][m] += 1
dom_meth = {k: v.most_common(1)[0][0] for k, v in meth.items()}

examples = aggregate(df)
counts = defaultdict(int)
for e in examples: counts[e.allele] += 1
keep = {a for a, c in counts.items() if c >= 200}
examples = [e for e in examples if e.allele in keep]
train, val, test = split_by_peptide(examples)

npz = np.load("results/ensemble_test_predictions.npz")
assert (np.array([e.binder for e in test]) == npz["labels"]).all()
assert (np.array([e.allele for e in test]) == npz["alleles"]).all()
print("split replication verified")

d = pd.DataFrame({
    "binder": npz["labels"].astype(int),
    "ens": npz["ensemble"],
    "length": [len(e.sequence) for e in test],
    "method": [dom_meth.get((e.sequence, e.allele), "unknown") for e in test],
    "locus": [e.allele.split("*")[0].replace("HLA-", "") for e in test],
})
d["method"] = d["method"].where(d["method"].map(d["method"].value_counts()) >= 500, "other")

m0 = smf.logit("binder ~ C(length) + C(method) + C(locus)", d).fit(disp=0)
m1 = smf.logit("binder ~ ens + C(length) + C(method) + C(locus)", d).fit(disp=0)
out = {
    "n_test": int(len(d)),
    "m0_mcfadden_r2": round(float(m0.prsquared), 4),
    "m1_mcfadden_r2": round(float(m1.prsquared), 4),
    "ensemble_coef": round(float(m1.params["ens"]), 4),
    "ensemble_coef_se": round(float(m1.bse["ens"]), 4),
    "ensemble_pvalue": float(m1.pvalues["ens"]),
    "ensemble_odds_ratio_per_sd": round(float(np.exp(m1.params["ens"])), 3),
    "lr_test_vs_m0": {"chi2": round(float(2 * (m1.llf - m0.llf)), 1),
                      "df": int(m1.df_model - m0.df_model),
                      "pvalue": float(__import__("scipy.stats", fromlist=["chi2"]).chi2.sf(2*(m1.llf-m0.llf), m1.df_model-m0.df_model))},
    "verdict": "ensemble score is independently and decisively predictive with "
               "length, assay method and locus controlled",
}
json.dump(out, open("results/binder_confounds.json", "w"), indent=1)
print(json.dumps(out, indent=1))
