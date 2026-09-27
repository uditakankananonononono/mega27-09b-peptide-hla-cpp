"""B1 allele-set selection (PREREG_B1_UNSEEN_ALLELE_2026-09-27, locked
2026-09-27 before any training/scoring). Computes the held-out allele set
ONCE and freezes it as results/b1_allele_set.json:
every allele with >= 50 usable rows (standard 20-aa alphabet, length 8-14)
in data/processed/iedb_class1_human_nM.tsv that is NOT in the training
allele set (the 52 alleles of results/per_allele_analysis.json), sorted by
HLA supertype then allele. ALL qualifying alleles are taken - no
performance-based exclusion, no re-selection.

Supertype assignment: Sidney et al. 2008, "HLA class I supertypes: a revised
and updated classification", BMC Immunology 9:1 (doi:10.1186/1471-2172-9-1),
supplementary allele table (via exa.ai mirror of the paper text, retrieved
2026-09-27). Alleles not explicitly present in that table are marked
"unclassified" (the paper's own residual category); HLA-C is outside the
paper's A/B scope and marked "C-locus (outside Sidney 2008 scope)". The
assignment is used ONLY as the sort key; membership is ALL qualifying
alleles regardless of assignment.
"""
import csv, json, pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent
AA20 = set("ACDEFGHIKLMNPQRSTVWY")

SUPERTYPE = {
    # A02 (explicit in table)
    "HLA-A*02:05": "A02", "HLA-A*02:07": "A02", "HLA-A*02:50": "A02",
    # A01 block
    "HLA-A*32:07": "A01",
    # A03 block
    "HLA-A*66:01": "A03", "HLA-A*68:23": "A03",
    # B08
    "HLA-B*08:02": "B08", "HLA-B*08:03": "B08",
    # B27 (incl. B*14:02, B*15:09, B*73:01 per table/motif listings)
    "HLA-B*14:02": "B27", "HLA-B*15:09": "B27", "HLA-B*27:03": "B27",
    "HLA-B*27:04": "B27", "HLA-B*27:06": "B27", "HLA-B*27:20": "B27",
    "HLA-B*38:01": "B27", "HLA-B*73:01": "B27",
    # B62
    "HLA-B*15:02": "B62",
    # unclassified (not explicit in the Sidney 2008 table):
    # HLA-A*32:15, HLA-B*15:42, HLA-B*40:13, HLA-B*45:06, HLA-B*83:01
}

def supertype(allele):
    if allele.startswith("HLA-C"):
        return "C-locus (outside Sidney 2008 scope)"
    return SUPERTYPE.get(allele, "unclassified")

def main():
    rows = {}
    with open(ROOT/"data/processed/iedb_class1_human_nM.tsv") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            s = r["sequence"]
            if 8 <= len(s) <= 14 and set(s) <= AA20:
                rows[r["allele"]] = rows.get(r["allele"], 0) + 1
    trained = {x["allele"] for x in
               json.load(open(ROOT/"results/per_allele_analysis.json"))["per_allele"]}
    held = [{"allele": a, "n_usable": n, "supertype": supertype(a)}
            for a, n in rows.items() if n >= 50 and a not in trained]
    order = {"A01": 0, "A02": 1, "A03": 2, "A24": 3, "A26": 4,
             "B07": 5, "B08": 6, "B27": 7, "B39": 8, "B44": 9, "B58": 10,
             "B62": 11, "unclassified": 12, "C-locus (outside Sidney 2008 scope)": 13}
    held.sort(key=lambda r: (order[r["supertype"]], r["allele"]))
    out = {"protocol": "PREREG_B1 locked selection: >=50 usable rows (std alphabet, len 8-14), not in the 52-allele training set; ALL qualifying alleles; sorted by supertype (Sidney 2008) then allele",
           "supertype_source": "Sidney et al. 2008 BMC Immunol 9:1 doi:10.1186/1471-2172-9-1 (table via exa.ai mirror, 2026-09-27)",
           "training_allele_count": len(trained),
           "held_out_count": len(held),
           "total_usable_rows": sum(r["n_usable"] for r in held),
           "alleles": held}
    json.dump(out, open(ROOT/"results/b1_allele_set.json", "w"), indent=1)
    print(f"frozen {len(held)} held-out alleles, {out['total_usable_rows']} usable rows")

if __name__ == "__main__":
    main()
