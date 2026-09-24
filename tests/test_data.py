"""Hermetic tests for data parsers and dataset construction."""
import math

import pytest

from peptidehlacpp.data.cppsite import (Peptide, dedup_exact, kmer_jaccard,
                                       natural_only, parse_fasta,
                                       redundancy_filter,
                                       sample_length_matched_windows)
from peptidehlacpp.data.iedb import normalize_allele, row_is_usable
from peptidehlacpp.data.iedb_dataset import PHLAExample, aggregate, split_by_peptide

FASTA = ">p1\nACDEFGHIK\n>p2\nklmnpqrst\n>p3\nACD*EFG\n"


def test_parse_fasta(tmp_path):
    p = tmp_path / "t.fa"
    p.write_text(FASTA)
    recs = parse_fasta(p)
    assert [r.identifier for r in recs] == ["p1", "p2", "p3"]
    assert recs[1].sequence == "KLMNPQRST"  # uppercased


def test_natural_and_dedup(tmp_path):
    p = tmp_path / "t.fa"
    p.write_text(FASTA + ">p4\nACDEFGHIK\n")
    recs = parse_fasta(p)
    nat = natural_only(recs)
    assert len(nat) == 3  # p3 has a stop char
    assert len(dedup_exact(nat)) == 2


def test_kmer_jaccard_bounds():
    assert kmer_jaccard("ACDEFG", "ACDEFG") == 1.0
    assert kmer_jaccard("AAAAAA", "CCCCCC") == 0.0
    assert 0 < kmer_jaccard("AAAAAC", "AAAAAG") < 1


def test_redundancy_filter_removes_near_dupes():
    recs = [Peptide("a", "ACDEFGHIK"), Peptide("b", "ACDEFGHIQ"),
            Peptide("c", "WYWVYWVYW")]
    kept = redundancy_filter(recs, threshold=0.6)
    names = {r.identifier for r in kept}
    assert "c" in names
    assert len(kept) == 2  # a and b collapse


def test_negative_windows_match_lengths():
    pos = [Peptide("p", "ACDEFGHIK"), Peptide("q", "ACDEFGHIKLMNPQRST")]
    pool = [Peptide(f"n{i}", "A" * 60 + "CDEFGH") for i in range(20)]
    neg = sample_length_matched_windows(pool, pos, n_per_pos=2, seed=1)
    assert len(neg) == 4
    assert sorted(len(n.sequence) for n in neg) == [9, 9, 17, 17]


def test_normalize_allele():
    assert normalize_allele("HLA-A*02:01") == "HLA-A*02:01"
    assert normalize_allele("human HLA-B*07:02 ref") == "HLA-B*07:02"
    assert normalize_allele("H-2-Kb") is None  # mouse, not HLA


def _fake_row(seq, allele, units, value, mhc_class="I"):
    row = [""] * 112
    row[11] = seq
    row[43] = "Homo sapiens"
    row[90] = "purified MHC/direct/fluorescence"
    row[92] = units
    row[94] = "Positive"
    row[95] = "="
    row[96] = str(value)
    row[107] = allele
    row[111] = mhc_class
    return row


def test_row_is_usable_filters():
    keep, seq, allele, units, val, qual, ref = row_is_usable(
        _fake_row("ACDEFGHIK", "HLA-A*02:01", "nM", "150.5"))
    assert keep and seq == "ACDEFGHIK" and allele == "HLA-A*02:01" and val == 150.5
    assert not row_is_usable(_fake_row("ACD", "HLA-A*02:01", "nM", "150"))[0]      # too short
    assert not row_is_usable(_fake_row("ACDEFGHIK", "HLA-A*02:01", "uM", "150"))[0]  # wrong units
    assert not row_is_usable(_fake_row("ACDEFGHIK", "HLA-A*02:01", "nM", "-5"))[0]   # bad value
    assert not row_is_usable(_fake_row("ACDEFGHIK", "HLA-A*02:01", "nM", "150", "II"))[0]  # class II
    assert not row_is_usable(_fake_row("ACDEXGHIK", "HLA-A*02:01", "nM", "150"))[0]  # non-AA


def test_aggregate_geometric_mean():
    import pandas as pd
    df = pd.DataFrame({
        "sequence": ["ACDEFGHIK", "ACDEFGHIK", "KLMNPQRST"],
        "allele": ["HLA-A*02:01"] * 3,
        "inequality": ["=", "=", "="],
        "value_nm": [10.0, 1000.0, 5.0],
    })
    ex = aggregate(df)
    assert len(ex) == 2
    first = next(e for e in ex if e.sequence == "ACDEFGHIK")
    assert math.isclose(first.ic50_nm, 100.0, rel_tol=1e-9)  # geo mean of 10, 1000
    assert first.binder is True
    assert first.n_measurements == 2


def test_split_has_no_peptide_leakage():
    rng_peps = [f"ACDEFGHI{chr(65 + i % 20)}" for i in range(200)]
    ex = [PHLAExample(sequence=p, allele="HLA-A*02:01", ic50_nm=100.0,
                      log_ic50=2.0, binder=True, n_measurements=1) for p in rng_peps]
    tr, va, te = split_by_peptide(ex, seed=1)
    s_tr = {e.sequence for e in tr}
    assert s_tr.isdisjoint({e.sequence for e in va})
    assert s_tr.isdisjoint({e.sequence for e in te})
    assert len(tr) + len(va) + len(te) == len(ex)
