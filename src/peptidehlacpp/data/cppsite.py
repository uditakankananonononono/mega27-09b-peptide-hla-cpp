"""Parsers for CPPsite 2.0 fasta files and UniProt negative sampling."""
from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

AA = set("ACDEFGHIKLMNPQRSTVWY")


@dataclass(frozen=True)
class Peptide:
    identifier: str
    sequence: str


def parse_fasta(path: str | Path) -> list[Peptide]:
    """Parse a fasta file into Peptide records, uppercased, whitespace stripped."""
    records: list[Peptide] = []
    name: str | None = None
    chunks: list[str] = []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append(Peptide(name, "".join(chunks).upper()))
                name = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.replace(" ", ""))
    if name is not None:
        records.append(Peptide(name, "".join(chunks).upper()))
    return records


def natural_only(records: list[Peptide]) -> list[Peptide]:
    """Keep sequences composed only of the 20 canonical amino acids."""
    return [r for r in records if r.sequence and set(r.sequence) <= AA]


def dedup_exact(records: list[Peptide]) -> list[Peptide]:
    seen: set[str] = set()
    out: list[Peptide] = []
    for r in records:
        if r.sequence not in seen:
            seen.add(r.sequence)
            out.append(r)
    return out


def kmer_jaccard(a: str, b: str, k: int = 3) -> float:
    """Jaccard similarity over k-mer sets - cheap homology proxy."""
    sa = {a[i:i + k] for i in range(len(a) - k + 1)} or {a}
    sb = {b[i:i + k] for i in range(len(b) - k + 1)} or {b}
    inter = len(sa & sb)
    return inter / (len(sa) + len(sb) - inter) if (sa or sb) else 0.0


def redundancy_filter(records: list[Peptide], threshold: float = 0.6,
                      k: int = 3) -> list[Peptide]:
    """Greedy redundancy removal: drop a record if its k-mer Jaccard similarity
    to any already-kept record exceeds `threshold`. O(n^2) worst case; fine for
    the ~1.5k CPPsite scale."""
    kept: list[Peptide] = []
    for r in records:
        if all(kmer_jaccard(r.sequence, q.sequence, k) < threshold for q in kept):
            kept.append(r)
    return kept


def sample_length_matched_windows(neg_pool: list[Peptide], positives: list[Peptide],
                                  n_per_pos: int = 1, seed: int = 13,
                                  max_len: int = 35) -> list[Peptide]:
    """Sample negative windows from longer UniProt sequences, matched to the
    positive length distribution (standard CPP-benchmark practice)."""
    rng = random.Random(seed)
    eligible = [p for p in neg_pool if len(p.sequence) >= 8]
    by_len: dict[int, list[Peptide]] = {}
    for p in eligible:
        by_len.setdefault(len(p), []).append(p)
    out: list[Peptide] = []
    for pos in positives:
        L = len(pos.sequence)
        for _ in range(n_per_pos):
            donors = by_len.get(L) or eligible
            donor = rng.choice(donors)
            if len(donor.sequence) == L:
                frag = donor.sequence
            else:
                start = rng.randrange(0, len(donor.sequence) - L + 1)
                frag = donor.sequence[start:start + L]
            if set(frag) <= AA and 8 <= len(frag) <= max_len:
                out.append(Peptide(f"neg|{donor.identifier}|{start if len(donor.sequence) != L else 0}", frag))
    return out
