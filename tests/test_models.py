"""Hermetic tests for model forward passes and the PSSM baseline."""
import numpy as np
import torch

from peptidehlacpp import features as F
from peptidehlacpp.models.cnn import PHLACNN
from peptidehlacpp.models.gnn import PHLAGNN
from peptidehlacpp.models.pssm import AllelePSSM


def test_cnn_forward_shapes():
    m = PHLACNN(n_alleles=4)
    x = torch.randn(3, 15, F.CNN_CHANNELS)
    mask = torch.ones(3, 15)
    al = torch.tensor([0, 1, 2])
    r, c = m(x, mask, al)
    assert r.shape == (3,) and c.shape == (3,)
    assert torch.isfinite(r).all() and torch.isfinite(c).all()


def test_gnn_forward_shapes_and_masking():
    m = PHLAGNN(n_alleles=4)
    x = torch.randn(2, 15, F.CNN_CHANNELS)
    mask = torch.zeros(2, 15); mask[:, :9] = 1.0
    adj = torch.from_numpy(np.stack([F.backbone_graph(9, 15)] * 2))
    al = torch.tensor([0, 3])
    r, c = m(x, adj, mask, al)
    assert r.shape == (2,)
    assert torch.isfinite(c).all()


def test_gnn_adj_normalization():
    adj = torch.from_numpy(F.backbone_graph(9, 15)).unsqueeze(0)
    a_hat = PHLAGNN.normalize_adj(adj)
    assert torch.isfinite(a_hat).all()
    assert (a_hat.diagonal(dim1=-2, dim2=-1) > 0).all()  # self loops


def test_pssm_recovers_additive_signal():
    rng = np.random.default_rng(0)
    W_true = rng.normal(0, 1, size=(15, 20))
    seqs = []
    ys = []
    for _ in range(400):
        s = "".join(rng.choice(list(F.AA_ORDER), size=9))
        y = sum(W_true[i, F.AA_INDEX[a]] for i, a in enumerate(s))
        seqs.append(s); ys.append(y)
    m = AllelePSSM(ridge=0.1).fit(seqs, np.array(ys))
    pred = m.predict(seqs[:50])
    corr = np.corrcoef(pred, ys[:50])[0, 1]
    assert corr > 0.95  # ridge fit recovers the additive ground truth
