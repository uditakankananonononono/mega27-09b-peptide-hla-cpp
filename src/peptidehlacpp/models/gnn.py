"""Allele-conditioned GNN over the peptide residue graph.

Message passing is implemented directly in torch (no torch_geometric):
GCN-style normalized aggregation A_hat = D^-1/2 (A + I) D^-1/2 with residual
updates. Node features = channel stack (44); graph = backbone + next-nearest
+ anchor edges (see features.backbone_graph). Readout: masked mean+max concat
with allele embedding -> MLP heads (log10 IC50 regression + binder logit).
"""
from __future__ import annotations

import torch
import torch.nn as nn

from ..features import CNN_CHANNELS, MAX_LEN_PHLA


class GCNLayer(nn.Module):
    def __init__(self, cin: int, cout: int):
        super().__init__()
        self.lin = nn.Linear(cin, cout)
        self.act = nn.GELU()

    def forward(self, h: torch.Tensor, a_hat: torch.Tensor,
                mask: torch.Tensor) -> torch.Tensor:
        out = self.act(torch.bmm(a_hat, self.lin(h)))
        return out * mask.unsqueeze(-1)


class PHLAGNN(nn.Module):
    def __init__(self, n_alleles: int, allele_dim: int = 32, hidden: int = 96,
                 layers: int = 3, max_len: int = MAX_LEN_PHLA, dropout: float = 0.15):
        super().__init__()
        self.allele_emb = nn.Embedding(n_alleles, allele_dim)
        self.inp = nn.Linear(CNN_CHANNELS, hidden)
        self.layers = nn.ModuleList(GCNLayer(hidden, hidden) for _ in range(layers))
        self.norm = nn.LayerNorm(hidden)
        pooled = 2 * hidden + allele_dim
        self.head = nn.Sequential(
            nn.Linear(pooled, 192), nn.GELU(), nn.Dropout(dropout),
            nn.Linear(192, 64), nn.GELU(),
        )
        self.reg_out = nn.Linear(64, 1)
        self.cls_out = nn.Linear(64, 1)

    @staticmethod
    def normalize_adj(adj: torch.Tensor) -> torch.Tensor:
        a = adj + torch.eye(adj.size(-1), device=adj.device)
        deg = a.sum(-1).clamp(min=1.0)
        d_inv_sqrt = deg.pow(-0.5)
        return a * d_inv_sqrt.unsqueeze(-1) * d_inv_sqrt.unsqueeze(-2)

    def forward(self, x: torch.Tensor, adj: torch.Tensor, mask: torch.Tensor,
                allele_idx: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        a_hat = self.normalize_adj(adj)
        h = self.inp(x) * mask.unsqueeze(-1)
        for layer in self.layers:
            h = self.norm(h + layer(h, a_hat, mask))
        m = mask.unsqueeze(-1)
        denom = m.sum(dim=1).clamp(min=1.0)
        mean_pool = (h * m).sum(1) / denom
        max_pool = h.masked_fill(m == 0, -1e4).amax(dim=1)
        pooled = torch.cat([mean_pool, max_pool, self.allele_emb(allele_idx)], dim=1)
        z = self.head(pooled)
        return self.reg_out(z).squeeze(-1), self.cls_out(z).squeeze(-1)
