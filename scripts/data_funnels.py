"""Dataset-construction funnel counts -> results/data_funnels.json."""
import json
import sys
sys.path.insert(0, "src")
import pandas as pd
from peptidehlacpp.data.cppsite import (parse_fasta, natural_only, dedup_exact,
                                        redundancy_filter)
from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv

raw = parse_fasta("data/raw/cppsite2_natural.fa")
nat = natural_only(raw)
dd = dedup_exact(nat)
rf = redundancy_filter(dd)
cpp_len = [p for p in rf if 8 <= len(p.sequence) <= 40]
df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
ex = aggregate(df)
alleles = {e.allele for e in ex}
from collections import Counter
cnt = Counter(e.allele for e in ex)
out = {
    "iedb": {"filtered_measurements": len(df), "unique_pairs": len(ex),
             "alleles": len(alleles),
             "alleles_ge200": sum(1 for a in alleles if cnt[a] >= 200)},
    "cppsite": {"raw": len(raw), "natural": len(nat), "dedup": len(dd),
                "redundancy_filtered": len(rf), "length_window": len(cpp_len)},
}
json.dump(out, open("results/data_funnels.json", "w"), indent=1)
print(json.dumps(out, indent=1))
