"""Allele-conditioned CNN for peptide-HLA binding (regression + classification).

Architecture: per-position channel stack (44 ch) -> 3 stacked 1D conv blocks
(kernels 3/5/7, dilated) with GELU + batch-norm -> masked mean+max pool ->
concat with allele embedding -> MLP head predicting log10(IC50); a sigmoid
head gives binder probability at 500 nM. Designed to fit in <1.5GB RAM.
"""
from __future__ import annotations

import torch
import torch.nn as nn

from ..features import CNN_CHANNELS, MAX_LEN_PHLA


class ConvBlock(nn.Module):
    def __init__(self, cin: int, cout: int, k: int, dilation: int = 1):
        super().__init__()
        pad = dilation * (k - 1) // 2
        self.conv = nn.Conv1d(cin, cout, k, padding=pad, dilation=dilation)
        self.bn = nn.BatchNorm1d(cout)
        self.act = nn.GELU()

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.act(self.bn(self.conv(x)))


class PHLACNN(nn.Module):
    def __init__(self, n_alleles: int, allele_dim: int = 32,
                 channels: tuple[int, ...] = (96, 96, 96),
                 max_len: int = MAX_LEN_PHLA, dropout: float = 0.15):
        super().__init__()
        self.max_len = max_len
        self.allele_emb = nn.Embedding(n_alleles, allele_dim)
        blocks = []
        cin = CNN_CHANNELS
        for i, cout in enumerate(channels):
            blocks.append(ConvBlock(cin, cout, k=3 + 2 * (i % 3), dilation=1 + i))
            cin = cout
        self.blocks = nn.ModuleList(blocks)
        pooled = 2 * channels[-1] + allele_dim
        self.head = nn.Sequential(
            nn.Linear(pooled, 192), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(192, 64), nn.GELU(),
        )
        self.reg_out = nn.Linear(64, 1)
        self.cls_out = nn.Linear(64, 1)

    def forward(self, x: torch.Tensor, mask: torch.Tensor,
                allele_idx: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        # x: (B, L, C) -> conv over L
        h = x.transpose(1, 2)
        for b in self.blocks:
            h = b(h)
        m = mask.unsqueeze(1)  # (B,1,L)
        h = h * m
        denom = m.sum(dim=2).clamp(min=1.0)
        mean_pool = h.sum(dim=2) / denom
        max_pool = h.masked_fill(m == 0, -1e4).amax(dim=2)
        pooled = torch.cat([mean_pool, max_pool, self.allele_emb(allele_idx)], dim=1)
        z = self.head(pooled)
        return self.reg_out(z).squeeze(-1), self.cls_out(z).squeeze(-1)
