"""Hermetic tests for the peptidehlacpp CLI (loads committed artifacts only)."""
import json
from pathlib import Path

import numpy as np
import pytest

from peptidehlacpp.cli import CPPScreener, PHLAPredictor, main, read_fasta

WDIR = Path("results")
pytestmark = pytest.mark.skipif(not (WDIR / "pssm_weights.npz").exists(),
                                reason="model artifacts not present")


def test_predict_known_epitope():
    pred = PHLAPredictor(WDIR)
    rows = pred.predict("HLA-A*02:01", ["GILGFVFTL", "AAAAAAAAA"])
    assert rows[0]["cnn_p_binder"] > 0.9  # flu MP epitope must rank as binder
    assert rows[1]["cnn_p_binder"] < 0.1  # poly-Ala must not
    assert rows[0]["pssm_ic50_nM"] < rows[1]["pssm_ic50_nM"]
    zs = [r["ensemble_z"] for r in rows]
    assert all(z is not None for z in zs) and abs(sum(zs)) < 1e-6


def test_predict_single_peptide_ensemble_nan(capsys):
    pred = PHLAPredictor(WDIR)
    rows = pred.predict("HLA-A*02:01", ["GILGFVFTL"])
    assert rows[0]["ensemble_z"] is None  # z-score undefined for batch of 1


def test_predict_allele_coverage_error():
    pred = PHLAPredictor(WDIR)
    with pytest.raises(ValueError, match="not covered"):
        pred.predict("HLA-Z*99:99", ["GILGFVFTL"])


def test_pssm_weights_match_refit():
    """Committed weights must reproduce a direct AllelePSSM computation."""
    from peptidehlacpp.models.pssm import AllelePSSM
    d = np.load(WDIR / "pssm_weights.npz")
    i = list(d["alleles"]).index("HLA-A*02:01")
    m = AllelePSSM()
    m.W, m.b = d["W"][i], float(d["b"][i])
    pred = PHLAPredictor(WDIR)
    seqs = ["SIINFEKL", "GILGFVFTL", "KLGEFYNQM"]
    assert np.allclose(pred.pssm_log10ic50("HLA-A*02:01", seqs), m.predict(seqs))


def test_screen_ranking_and_determinism():
    scr = CPPScreener(WDIR)
    seqs = ["RQIKIWFQNRRMKWKK", "AAAAAAAAAAAA"]
    r1 = scr.screen(seqs)
    r2 = scr.screen(seqs)
    assert r1 == r2  # deterministic
    by_seq = {r["sequence"]: r for r in r1}
    assert by_seq["RQIKIWFQNRRMKWKK"]["p_cpp"] > 0.9  # penetratin
    assert by_seq["AAAAAAAAAAAA"]["p_cpp"] < 0.5


def test_design_cpp_seeded_valid(tmp_path):
    rc = main(["--weights-dir", str(WDIR), "design-cpp", "--n", "3",
               "--sample", "300", "--seed", "5", "--json"])
    assert rc == 0


def test_design_cpp_output_properties(capsys):
    main(["--weights-dir", str(WDIR), "design-cpp", "--n", "3",
          "--sample", "300", "--seed", "5", "--json"])
    rows = json.loads(capsys.readouterr().out)
    train = set(Path("results/cpp_train_sequences.txt").read_text().split())
    assert 1 <= len(rows) <= 3
    for r in rows:
        assert 8 <= len(r["sequence"]) <= 40
        assert r["sequence"] not in train  # novelty screen enforced
        assert r["p_cpp"] >= 0.7           # cascade classifier gate
        assert 2.0 <= r["net_charge"] <= 12.0


def test_read_fasta(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(">a desc\nACD\nEF\n>b\nGG\n")
    assert read_fasta(str(f)) == [("a", "ACDEF"), ("b", "GG")]


def test_cli_predict_json_e2e(capsys):
    rc = main(["--weights-dir", str(WDIR), "predict", "--allele", "HLA-A*02:01",
               "--peptides", "GILGFVFTL", "AAAAAAAAA", "--json"])
    assert rc == 0
    rows = json.loads(capsys.readouterr().out)
    assert len(rows) == 2 and rows[0]["binder_500nM"] and not rows[1]["binder_500nM"]
