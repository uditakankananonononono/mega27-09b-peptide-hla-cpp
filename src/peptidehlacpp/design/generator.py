"""Autoregressive GRU peptide generator trained on CPPsite CPPs, with an
in-silico screening cascade: classifier score, net charge, hydrophobic
moment, mean hydropathy, and sequence-novelty check against the training set.

The generator learns p(x_t | x_<t) over amino-acid tokens; sampling is
ancestral with temperature. Screens are real biophysical filters used in the
CPP literature (cationicity + amphipathicity window).
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn

from ..features import (AA_ORDER, hydrophobic_moment, mean_hydropathy,
                        net_charge)

VOCAB = ["<s>", "<e>"] + list(AA_ORDER)
TOK = {t: i for i, t in enumerate(VOCAB)}
MAX_GEN_LEN = 40


class CPPGenerator(nn.Module):
    def __init__(self, emb: int = 32, hidden: int = 128, layers: int = 1):
        super().__init__()
        self.emb = nn.Embedding(len(VOCAB), emb)
        self.gru = nn.GRU(emb, hidden, num_layers=layers, batch_first=True)
        self.out = nn.Linear(hidden, len(VOCAB))

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        h, _ = self.gru(self.emb(tokens))
        return self.out(h)


def encode_seq(seq: str, max_len: int = MAX_GEN_LEN) -> list[int]:
    toks = [TOK["<s>"]] + [TOK[a] for a in seq[:max_len]] + [TOK["<e>"]]
    return toks


def train_generator(seqs: list[str], epochs: int = 30, lr: float = 2e-3,
                    batch: int = 128, seed: int = 9, log_every: int = 5) -> CPPGenerator:
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    model = CPPGenerator()
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-5)
    ce = nn.CrossEntropyLoss(ignore_index=-100)
    enc = [encode_seq(s) for s in seqs]
    n = len(enc)
    for ep in range(epochs):
        model.train()
        order = rng.permutation(n)
        total, count = 0.0, 0
        for s in range(0, n, batch):
            chunk = [enc[i] for i in order[s:s + batch]]
            L = max(len(c) for c in chunk)
            inp = torch.full((len(chunk), L - 1), TOK["<e>"], dtype=torch.long)
            tgt = torch.full((len(chunk), L - 1), -100, dtype=torch.long)
            for r, c in enumerate(chunk):
                inp[r, :len(c) - 1] = torch.tensor(c[:-1])
                tgt[r, :len(c) - 1] = torch.tensor(c[1:])
            logits = model(inp)
            loss = ce(logits.reshape(-1, len(VOCAB)), tgt.reshape(-1))
            opt.zero_grad(); loss.backward(); opt.step()
            total += loss.item() * len(chunk); count += len(chunk)
        if (ep + 1) % log_every == 0:
            ppl = float(np.exp(total / max(1, count)))
            print(f"gen epoch {ep+1}/{epochs} perplexity={ppl:.2f}", flush=True)
    return model


@torch.no_grad()
def sample(model: CPPGenerator, n: int, temperature: float = 0.9,
           seed: int = 17) -> list[str]:
    model.eval()
    g = torch.Generator().manual_seed(seed)
    out: list[str] = []
    while len(out) < n:
        tokens = torch.tensor([[TOK["<s>"]]], dtype=torch.long)
        chars: list[str] = []
        for _ in range(MAX_GEN_LEN):
            logits = model(tokens)[0, -1] / temperature
            probs = torch.softmax(logits, dim=-1)
            nxt = torch.multinomial(probs, 1, generator=g).item()
            sym = VOCAB[nxt]
            if sym == "<e>":
                break
            if sym == "<s>":
                continue
            chars.append(sym)
            tokens = torch.cat([tokens, torch.tensor([[nxt]])], dim=1)
        seq = "".join(chars)
        if 8 <= len(seq) <= MAX_GEN_LEN:
            out.append(seq)
    return out


def screen_candidates(candidates: list[str], classifier_score,
                      train_seqs: set[str],
                      min_score: float = 0.7,
                      charge_range: tuple[float, float] = (2.0, 12.0),
                      min_moment: float = 0.15,
                      max_hydropathy: float = 1.5) -> list[dict]:
    """Apply the in-silico cascade. classifier_score: seq -> p(CPP)."""
    passed: list[dict] = []
    seen: set[str] = set()
    for s in candidates:
        if s in seen or s in train_seqs:
            continue
        seen.add(s)
        q = net_charge(s)
        mu = hydrophobic_moment(s)
        h = mean_hydropathy(s)
        if not (charge_range[0] <= q <= charge_range[1]):
            continue
        if mu < min_moment or h > max_hydropathy:
            continue
        p = classifier_score(s)
        if p < min_score:
            continue
        passed.append({"sequence": s, "p_cpp": round(float(p), 4),
                       "net_charge": round(q, 2),
                       "hydrophobic_moment": round(mu, 3),
                       "mean_hydropathy": round(h, 3)})
    passed.sort(key=lambda d: -d["p_cpp"])
    return passed
