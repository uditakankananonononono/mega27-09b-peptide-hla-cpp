"""Peptide encodings: one-hot, BLOSUM62, physicochemical scales, graph views."""
from __future__ import annotations

import numpy as np

AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
AA_INDEX = {a: i for i, a in enumerate(AA_ORDER)}
PAD_INDEX = 20

# BLOSUM62 in AA_ORDER (from NCBI blastp distribution, public domain matrix)
_BLOSUM62_RAW = """
   A  C  D  E  F  G  H  I  K  L  M  N  P  Q  R  S  T  V  W  Y
A  4  0 -2 -1 -2  0 -2 -1 -1 -1 -1 -2 -1 -1 -1  1  0  0 -3 -2
C  0  9 -3 -4 -2 -3 -3 -1 -3 -1 -1 -3 -3 -3 -3 -1 -1 -1 -2 -2
D -2 -3  6  2 -3 -1 -1 -3 -1 -4 -3  1 -1  0 -2  0 -1 -3 -4 -3
E -1 -4  2  5 -3 -2  0 -3  1 -3 -2  0 -1  2  0  0 -1 -2 -3 -2
F -2 -2 -3 -3  6 -3 -1  0 -3  0  0 -3 -4 -3 -3 -2 -2 -1  1  3
G  0 -3 -1 -2 -3  6 -2 -4 -2 -4 -3  0 -2 -2 -2  0 -2 -3 -2 -3
H -2 -3 -1  0 -1 -2  8 -3 -1 -3 -2  1 -2  0  0 -1 -2 -3 -2  2
I -1 -1 -3 -3  0 -4 -3  4 -3  2  1 -3 -3 -3 -3 -2 -1  3 -3 -1
K -1 -3 -1  1 -3 -2 -1 -3  5 -2 -1  0 -1  1  2  0 -1 -2 -3 -2
L -1 -1 -4 -3  0 -4 -3  2 -2  4  2 -3 -3 -2 -2 -2 -1  1 -2 -1
M -1 -1 -3 -2  0 -3 -2  1 -1  2  5 -2 -2  0 -1 -1 -1  1 -1 -1
N -2 -3  1  0 -3  0  1 -3  0 -3 -2  6 -2  0  0  1  0 -3 -4 -2
P -1 -3 -1 -1 -4 -2 -2 -3 -1 -3 -2 -2  7 -1 -2 -1 -1 -2 -4 -3
Q -1 -3  0  2 -3 -2  0 -3  1 -2  0  0 -1  5  1  0 -1 -2 -2 -1
R -1 -3 -2  0 -3 -2  0 -3  2 -2 -1  0 -2  1  5 -1 -1 -3 -3 -2
S  1 -1  0  0 -2  0 -1 -2  0 -2 -1  1 -1  0 -1  4  1 -2 -3 -2
T  0 -1 -1 -1 -2 -2 -2 -2 -1 -1 -1 -1 -1 -2 -1 -1  5 -2 -2  0
V  0 -1 -3 -2 -1 -3 -2  3 -2  1  1 -3 -2 -2 -3 -2  0  4 -3 -1
W -3 -2 -4 -3  1 -2 -2 -3 -3 -2 -1 -4 -4 -2 -3 -3 -2 -3 11  2
Y -2 -2 -3 -2  3 -3  2 -1 -2 -1 -1 -2 -3 -1 -2 -2 -2 -1  2  7
"""

def _parse_blosum() -> np.ndarray:
    lines = [l.split() for l in _BLOSUM62_RAW.strip().splitlines()]
    header = lines[0]
    assert header == list(AA_ORDER)
    mat = np.zeros((20, 20), dtype=np.float32)
    for row in lines[1:]:
        i = AA_INDEX[row[0]]
        mat[i] = [float(x) for x in row[1:]]
    return mat

BLOSUM62 = _parse_blosum()

# Physicochemical scales (Kyte-Doolittle hydropathy; bulkiness (Zimmerman);
# net charge at pH 7; helix propensity (Pace-Scholtz)) - published constants.
HYDROPATHY = dict(zip(AA_ORDER,
    [1.8, 2.5, -3.5, -3.5, 2.8, -0.4, -3.2, 4.5, -3.9, 3.8,
     1.9, -3.5, -1.6, -3.5, -4.5, -0.8, -0.7, 4.2, -0.9, -1.3]))
CHARGE = dict(zip(AA_ORDER,
    [0, 0, -1, -1, 0, 0, 0.1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]))
HELIX_PROP = dict(zip(AA_ORDER,
    [1.42, 0.70, 1.01, 1.51, 1.13, 0.57, 1.00, 1.08, 1.16, 1.21,
     1.45, 0.67, 0.57, 1.11, 0.98, 0.77, 0.83, 1.06, 1.08, 0.69]))

MAX_LEN_PHLA = 15
MAX_LEN_CPP = 61


def one_hot(seq: str, max_len: int) -> np.ndarray:
    """(max_len, 21) float32 one-hot with a pad channel."""
    out = np.zeros((max_len, 21), dtype=np.float32)
    for i, a in enumerate(seq[:max_len]):
        out[i, AA_INDEX[a]] = 1.0
    for i in range(len(seq), max_len):
        out[i, PAD_INDEX] = 1.0
    return out


def blosum_enc(seq: str, max_len: int) -> np.ndarray:
    """(max_len, 20) BLOSUM62 rows, zeros at pads."""
    out = np.zeros((max_len, 20), dtype=np.float32)
    for i, a in enumerate(seq[:max_len]):
        out[i] = BLOSUM62[AA_INDEX[a]]
    return out


def physico_enc(seq: str, max_len: int) -> np.ndarray:
    """(max_len, 3) hydropathy/charge/helix-propensity, zeros at pads."""
    out = np.zeros((max_len, 3), dtype=np.float32)
    for i, a in enumerate(seq[:max_len]):
        out[i] = [HYDROPATHY[a], CHARGE[a], HELIX_PROP[a]]
    return out


def stacked_enc(seq: str, max_len: int) -> np.ndarray:
    """(max_len, 21+20+3) channel stack for CNN input."""
    return np.concatenate([one_hot(seq, max_len), blosum_enc(seq, max_len),
                           physico_enc(seq, max_len)], axis=1)


CNN_CHANNELS = 21 + 20 + 3


def backbone_graph(seq_len: int, max_len: int) -> np.ndarray:
    """Adjacency (max_len, max_len) for the residue graph: backbone edges
    (i,i+1) plus next-nearest (i,i+2) plus MHC-I anchor edges (2 <-> C-term),
    symmetrized. Self-loops added by the model's GCN normalization."""
    adj = np.zeros((max_len, max_len), dtype=np.float32)
    L = min(seq_len, max_len)
    for i in range(L - 1):
        adj[i, i + 1] = adj[i + 1, i] = 1.0
    for i in range(L - 2):
        adj[i, i + 2] = adj[i + 2, i] = 1.0
    if L >= 3:
        adj[1, L - 1] = adj[L - 1, 1] = 1.0  # P2 <-> P-Omega anchor coupling
    return adj


def hydrophobic_moment(seq: str, angle_deg: float = 100.0) -> float:
    """Eisenberg consensus hydrophobic moment for an alpha-helical stretch."""
    if not seq:
        return 0.0
    angle = np.deg2rad(angle_deg)
    vx = sum(HYDROPATHY[a] * np.cos(i * angle) for i, a in enumerate(seq))
    vy = sum(HYDROPATHY[a] * np.sin(i * angle) for i, a in enumerate(seq))
    return float(np.hypot(vx, vy) / len(seq))


def net_charge(seq: str) -> float:
    return float(sum(CHARGE[a] for a in seq))


def mean_hydropathy(seq: str) -> float:
    return float(np.mean([HYDROPATHY[a] for a in seq])) if seq else 0.0
