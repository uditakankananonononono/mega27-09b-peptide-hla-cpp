"""Accession-level dataset manifest: one entry per IEDB reference (study)
contributing to the retained class-I binding dataset, per the program's
dataset-counting rule (separate studies count separately; one study's
condition matrix is one dataset). Saves results/iedb_reference_manifest.json.
"""
import json
import pandas as pd

df = pd.read_csv("data/processed/iedb_class1_human_nM.tsv", sep="\t")
g = df.groupby("reference_iri").agg(
    n_measurements=("sequence", "size"),
    n_unique_peptides=("sequence", "nunique"),
    n_alleles=("allele", "nunique")).reset_index()
g = g.sort_values("n_measurements", ascending=False)
out = {
    "n_reference_datasets": int(len(g)),
    "total_measurements": int(g["n_measurements"].sum()),
    "top10": g.head(10).to_dict("records"),
    "references": g.to_dict("records"),
}
json.dump(out, open("results/iedb_reference_manifest.json", "w"), indent=1)
print(out["n_reference_datasets"], "references,",
      out["total_measurements"], "measurements; top:",
      [(r["reference_iri"].rsplit("/",1)[-1], r["n_measurements"]) for r in out["top10"][:5]])
