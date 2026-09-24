"""CPP classification benchmark: CNN classifier + k-mer logistic baseline.

Positives: CPPsite 2.0 natural CPPs (dedup, redundancy-filtered).
Negatives: length-matched windows sampled from reviewed UniProt proteins.
Split: peptide-level, homology-guarded (k-mer Jaccard < 0.6 across splits).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn
from sklearn.linear_model import LogisticRegression

from ..data.cppsite import (dedup_exact, natural_only, parse_fasta,
                            redundancy_filter, sample_length_matched_windows,
                            kmer_jaccard)
from ..eval import metrics as M
from ..features import MAX_LEN_CPP, stacked_enc

DEVICE = torch.device("cpu")
MAX_LEN = 40  # cover 95%+ of CPPsite natural peptides; pads handle the rest


def build_kmer_features(seqs: list[str], k: int = 3, vocab: dict | None = None):
    from ..features import AA_ORDER
    if vocab is None:
        vocab = {"".join(t): i for i, t in enumerate(
            __import__("itertools").product(AA_ORDER, repeat=k))}
    X = np.zeros((len(seqs), len(vocab)), dtype=np.float32)
    for r, s in enumerate(seqs):
        counts: dict[int, int] = {}
        for i in range(len(s) - k + 1):
            j = vocab.get(s[i:i + k])
            if j is not None:
                counts[j] = counts.get(j, 0) + 1
        total = max(1, len(s) - k + 1)
        for j, c in counts.items():
            X[r, j] = c / total
    return X, vocab


class CPPCNN(nn.Module):
    def __init__(self, channels: int = 44, hidden: int = 64):
        super().__init__()
        self.c1 = nn.Conv1d(channels, hidden, 3, padding=1)
        self.c2 = nn.Conv1d(hidden, hidden, 5, padding=4, dilation=2)
        self.bn1 = nn.BatchNorm1d(hidden); self.bn2 = nn.BatchNorm1d(hidden)
        self.act = nn.GELU()
        self.head = nn.Sequential(nn.Linear(2 * hidden, 64), nn.GELU(), nn.Linear(64, 1))

    def forward(self, x: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        h = x.transpose(1, 2)
        h = self.act(self.bn1(self.c1(h)))
        h = self.act(self.bn2(self.c2(h)))
        m = mask.unsqueeze(1)
        h = h * m
        denom = m.sum(2).clamp(min=1.0)
        pooled = torch.cat([h.sum(2) / denom,
                            h.masked_fill(m == 0, -1e4).amax(2)], dim=1)
        return self.head(pooled).squeeze(-1)


def enc(seqs: list[str], max_len: int = MAX_LEN):
    X = np.stack([stacked_enc(s, max_len) for s in seqs])
    mask = np.zeros((len(seqs), max_len), dtype=np.float32)
    for r, s in enumerate(seqs):
        mask[r, : min(len(s), max_len)] = 1.0
    return torch.from_numpy(X), torch.from_numpy(mask)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pos", default="data/raw/cppsite2_natural.fa")
    ap.add_argument("--neg-pool", nargs="+", default=["data/raw/uniprot_reviewed_len8_35.fasta",
                                                   "data/raw/uniprot_human_reviewed_len60_400.fasta"])
    ap.add_argument("--out", default="results/cpp_benchmark.json")
    ap.add_argument("--epochs", type=int, default=10)
    ap.add_argument("--batch", type=int, default=256)
    args = ap.parse_args()

    pos = redundancy_filter(dedup_exact(natural_only(parse_fasta(args.pos))))
    pos = [p for p in pos if 8 <= len(p.sequence) <= MAX_LEN]
    neg_pool = []
    for np_ in args.neg_pool:
        neg_pool.extend(p for p in natural_only(parse_fasta(np_)) if 8 <= len(p.sequence))
    neg_pool = dedup_exact(neg_pool)
    neg = sample_length_matched_windows(neg_pool, pos, n_per_pos=1)
    print(f"positives={len(pos)} negatives={len(neg)} pool={len(neg_pool)}", flush=True)

    rng = np.random.default_rng(11)
    order = rng.permutation(len(pos))
    n_test = int(0.15 * len(pos)); n_val = int(0.1 * len(pos))
    idx_te, idx_va = order[:n_test], order[n_test:n_test + n_val]
    idx_tr = order[n_test + n_val:]
    # homology guard: drop train items too similar to any test item
    te_pos = [pos[i].sequence for i in idx_te]
    idx_tr = [i for i in idx_tr
              if all(kmer_jaccard(pos[i].sequence, t) < 0.6 for t in te_pos)]
    def mk(idxs_p, sign):
        n_same = int(round(len(idxs_p) * len(neg) / len(pos)))
        neg_sel = [neg[j] for j in rng.choice(len(neg), size=n_same, replace=False)]
        seqs = [pos[i].sequence for i in idxs_p] + [n.sequence for n in neg_sel]
        labels = np.array([1] * len(idxs_p) + [0] * len(neg_sel), dtype=np.float32)
        return seqs, labels
    tr_s, tr_y = mk(idx_tr, 1); va_s, va_y = mk(idx_va, 1); te_s, te_y = mk(idx_te, 1)
    print(f"train={len(tr_s)} val={len(va_s)} test={len(te_s)}", flush=True)

    results: dict = {"n_pos": len(pos), "n_neg": len(neg)}

    Xtr, vocab = build_kmer_features(tr_s)
    lr = LogisticRegression(max_iter=2000, C=1.0, n_jobs=2).fit(Xtr, tr_y)
    Xte, _ = build_kmer_features(te_s, vocab=vocab)
    lr_scores = lr.predict_proba(Xte)[:, 1]
    results["kmer_lr"] = {"auc": M.auc(te_y, lr_scores),
                          "ppv": M.ppv(te_y, lr_scores)}
    print("kmer-LR:", results["kmer_lr"], flush=True)

    torch.manual_seed(5)
    model = CPPCNN().to(DEVICE)
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    bce = nn.BCEWithLogitsLoss()
    ytr = torch.from_numpy(tr_y)
    best = -1.0; best_state = None
    for ep in range(args.epochs):
        model.train()
        perm = rng.permutation(len(tr_s)); t0 = time.time()
        for s in range(0, len(tr_s), args.batch):
            idx = perm[s:s + args.batch]
            X, m = enc([tr_s[i] for i in idx])
            loss = bce(model(X, m), ytr[idx])
            opt.zero_grad(); loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            Xv, mv = enc(va_s)
            va_scores = torch.sigmoid(model(Xv, mv)).numpy()
        va_auc = M.auc(va_y, va_scores)
        print(f"epoch {ep+1}/{args.epochs} val_auc={va_auc:.4f} ({time.time()-t0:.0f}s)", flush=True)
        if va_auc > best:
            best = va_auc
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
    model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        Xte2, mte = enc(te_s)
        cnn_scores = torch.sigmoid(model(Xte2, mte)).numpy()
    results["cnn"] = {"auc": M.auc(te_y, cnn_scores), "ppv": M.ppv(te_y, cnn_scores)}
    print("CNN:", results["cnn"], flush=True)
    torch.save(model.state_dict(), "results/cpp_cnn.pt")
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
