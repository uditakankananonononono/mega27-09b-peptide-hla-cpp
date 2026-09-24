"""Turn the filtered IEDB TSV into a modeling dataset.

Aggregation rule (matches MHCflurry/NetMHCpan practice): for replicated
(peptide, allele) measurements, take the geometric mean of IC50 values and
the median inequality direction; qualitative labels only fill in where no
quantitative value exists for that pair.
"""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass

import numpy as np
import pandas as pd

BINDER_THRESHOLD_NM = 500.0


@dataclass(frozen=True)
class PHLAExample:
    sequence: str
    allele: str
    ic50_nm: float
    log_ic50: float
    binder: bool
    n_measurements: int


def load_filtered_tsv(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    df = df[df["inequality"].isin(["", "=", "<", ">", "<=", ">="]) | df["inequality"].isna()]
    df = df[df["value_nm"] > 0]
    return df


def aggregate(df: pd.DataFrame) -> list[PHLAExample]:
    """Geometric-mean aggregation of replicate (sequence, allele) assays."""
    groups: dict[tuple[str, str], list[float]] = defaultdict(list)
    for seq, allele, val in df[["sequence", "allele", "value_nm"]].itertuples(index=False):
        groups[(seq, allele)].append(float(val))
    examples: list[PHLAExample] = []
    for (seq, allele), vals in groups.items():
        geo = math.exp(sum(math.log(v) for v in vals) / len(vals))
        examples.append(PHLAExample(
            sequence=seq, allele=allele, ic50_nm=geo,
            log_ic50=math.log10(geo), binder=geo <= BINDER_THRESHOLD_NM,
            n_measurements=len(vals)))
    return examples


def allele_table(examples: list[PHLAExample]) -> pd.DataFrame:
    counts: dict[str, int] = defaultdict(int)
    for e in examples:
        counts[e.allele] += 1
    return (pd.DataFrame(sorted(counts.items(), key=lambda kv: -kv[1]),
                         columns=["allele", "n_measurements"]))


def split_by_peptide(examples: list[PHLAExample], val_frac: float = 0.1,
                     test_frac: float = 0.15, seed: int = 7
                     ) -> tuple[list[PHLAExample], list[PHLAExample], list[PHLAExample]]:
    """Peptide-level split: no peptide appears in more than one split (prevents
    the leakage that inflates random splits on this data)."""
    rng = np.random.default_rng(seed)
    pep_to_idx: dict[str, list[int]] = defaultdict(list)
    for i, e in enumerate(examples):
        pep_to_idx[e.sequence].append(i)
    peps = np.array(sorted(pep_to_idx))
    rng.shuffle(peps)
    n = len(peps)
    n_test = int(round(test_frac * n))
    n_val = int(round(val_frac * n))
    test_peps = set(peps[:n_test])
    val_peps = set(peps[n_test:n_test + n_val])
    train, val, test = [], [], []
    for pep, idxs in pep_to_idx.items():
        bucket = test if pep in test_peps else val if pep in val_peps else train
        bucket.extend(examples[i] for i in idxs)
    return train, val, test
