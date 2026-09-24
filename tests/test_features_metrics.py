"""Hermetic tests for encodings, graphs, and metrics."""
import numpy as np

from peptidehlacpp import features as F
from peptidehlacpp.eval import metrics as M


def test_one_hot_shape_and_pad():
    oh = F.one_hot("AC", 5)
    assert oh.shape == (5, 21)
    assert oh[0, F.AA_INDEX["A"]] == 1
    assert oh[4, F.PAD_INDEX] == 1


def test_blosum_diagonal_positive():
    enc = F.blosum_enc("W", 3)
    assert enc[0, F.AA_INDEX["W"]] == 11.0  # BLOSUM62 W-W score
    assert enc[2].sum() == 0.0  # pad row


def test_physico_known_values():
    enc = F.physico_enc("IK", 2)
    assert np.isclose(enc[0, 0], 4.5)  # Ile hydropathy
    assert np.isclose(enc[1, 2], 1.16)  # Lys helix propensity


def test_graph_anchor_edge():
    adj = F.backbone_graph(9, 15)
    assert adj[1, 8] == 1.0 and adj[8, 1] == 1.0  # P2 <-> C-term anchor
    assert adj[0, 1] == 1.0  # backbone
    assert (adj == adj.T).all()  # symmetric
    assert adj[10:, :].sum() == 0  # pad region empty


def test_hydrophobic_moment_against_formula():
    # Hand-computed: mu = |H_A + H_C*e^{i*100deg}| / 2
    angle = np.deg2rad(100.0)
    expect = abs(F.HYDROPATHY["A"] + F.HYDROPATHY["C"] * np.exp(1j * angle)) / 2
    assert np.isclose(F.hydrophobic_moment("AC"), expect)
    # hydrophobics aligned at the 100-degree period -> larger moment than scrambled
    aligned = "IDSSIDSSIDSS"
    scrambled = "IIDDDSSSSIIS"
    assert F.hydrophobic_moment(aligned) > F.hydrophobic_moment(scrambled)
    assert F.hydrophobic_moment("") == 0.0


def test_metrics_perfect_and_worst():
    labels = np.array([1, 1, 0, 0])
    good = np.array([0.9, 0.8, 0.2, 0.1])
    bad = np.array([0.1, 0.2, 0.8, 0.9])
    assert M.auc(labels, good) == 1.0
    assert M.auc(labels, bad) == 0.0
    assert M.ppv(labels, good, top_n=2) == 1.0
    assert M.srcc(np.array([1.0, 2.0, 3.0]), np.array([1.0, 2.0, 3.0])) == 1.0
    assert M.rmse(np.array([1.0]), np.array([1.0])) == 0.0


def test_auc_top_frac():
    labels = np.array([1] * 10 + [0] * 90)
    scores = np.concatenate([np.linspace(1, 0.5, 10), np.linspace(0.4, 0.1, 90)])
    assert 0.99 <= M.auc_top_frac(labels, scores, frac=0.1) <= 1.0  # all pos on top
    scores_bad = np.concatenate([np.linspace(0.1, 0.4, 10), np.linspace(0.5, 1, 90)])
    assert M.auc_top_frac(labels, scores_bad, frac=0.1) < 0.2
