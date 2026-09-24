"""Train PSSM / CNN / GNN on the filtered IEDB class-I dataset and evaluate.

Pipeline: aggregate replicates -> restrict to alleles with >= MIN_PER_ALLELE
measurements -> peptide-level split -> fit PSSM per allele; train CNN and GNN
pooled across alleles (allele-conditioned) -> metrics per model:
AUC(500nM), AUC0.1, PPV, SRCC on the held-out test split.
"""
from __future__ import annotations

import argparse
import json
import math
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import torch
import torch.nn as nn

from ..data.iedb_dataset import (BINDER_THRESHOLD_NM, PHLAExample, aggregate,
                                 load_filtered_tsv, split_by_peptide)
from ..eval import metrics as M
from ..features import MAX_LEN_PHLA, backbone_graph, stacked_enc
from ..models.cnn import PHLACNN
from ..models.gnn import PHLAGNN
from ..models.pssm import AllelePSSM

MIN_PER_ALLELE = 200
DEVICE = torch.device("cpu")


def enc_batch(seqs: list[str], allele_idx: np.ndarray, allele_map: dict[str, int],
              with_graph: bool):
    X = np.stack([stacked_enc(s, MAX_LEN_PHLA) for s in seqs])
    mask = np.zeros((len(seqs), MAX_LEN_PHLA), dtype=np.float32)
    for r, s in enumerate(seqs):
        mask[r, : min(len(s), MAX_LEN_PHLA)] = 1.0
    al = np.array([allele_map[a] for a in allele_idx], dtype=np.int64)
    X_t, mask_t, al_t = (torch.from_numpy(X), torch.from_numpy(mask),
                         torch.from_numpy(al))
    if with_graph:
        adj = np.stack([backbone_graph(len(s), MAX_LEN_PHLA) for s in seqs])
        return X_t, torch.from_numpy(adj), mask_t, al_t  # GNN forward order
    return X_t, mask_t, al_t  # CNN forward order


def train_torch_model(model, train: list[PHLAExample], val: list[PHLAExample],
                      allele_map: dict[str, int], epochs: int, batch: int,
                      lr: float, with_graph: bool, seed: int = 3):
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    opt = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    bce = nn.BCEWithLogitsLoss()
    mse = nn.MSELoss()
    y_reg = np.array([e.log_ic50 for e in train], dtype=np.float32)
    y_cls = np.array([e.binder for e in train], dtype=np.float32)
    y_reg_t, y_cls_t = torch.from_numpy(y_reg), torch.from_numpy(y_cls)
    n = len(train)
    best_val = math.inf
    best_state = None
    for ep in range(epochs):
        model.train()
        order = rng.permutation(n)
        t0 = time.time()
        for s in range(0, n, batch):
            idx = order[s:s + batch]
            seqs = [train[i].sequence for i in idx]
            als = np.array([train[i].allele for i in idx])
            tensors = enc_batch(seqs, als, allele_map, with_graph)
            pred_r, pred_c = model(*tensors)
            loss = mse(pred_r, y_reg_t[idx]) + bce(pred_c, y_cls_t[idx])
            opt.zero_grad(); loss.backward(); opt.step()
        # quick val AUC check
        model.eval()
        with torch.no_grad():
            va = evaluate_torch(model, val, allele_map, with_graph, batch=2048)
        print(f"epoch {ep+1}/{epochs} val_auc={va['auc']:.4f} "
              f"val_srcc={va['srcc']:.4f} ({time.time()-t0:.0f}s)", flush=True)
        if -va["auc"] < best_val:
            best_val = -va["auc"]
            best_state = {k: v.clone() for k, v in model.state_dict().items()}
    if best_state:
        model.load_state_dict(best_state)
    return model


@torch.no_grad()
def evaluate_torch(model, examples: list[PHLAExample], allele_map: dict[str, int],
                   with_graph: bool, batch: int = 2048) -> dict:
    model.eval()
    scores, regs = [], []
    for s in range(0, len(examples), batch):
        chunk = examples[s:s + batch]
        seqs = [e.sequence for e in chunk]
        als = np.array([e.allele for e in chunk])
        tensors = enc_batch(seqs, als, allele_map, with_graph)
        pr, pc = model(*tensors)
        regs.append(pr.numpy()); scores.append(pc.numpy())
    scores = np.concatenate(scores); regs = np.concatenate(regs)
    labels = np.array([e.binder for e in examples])
    true_reg = np.array([e.log_ic50 for e in examples])
    return {"auc": M.auc(labels, scores), "auc0.1": M.auc_top_frac(labels, scores),
            "ppv": M.ppv(labels, scores), "srcc": M.srcc(true_reg, regs),
            "rmse": M.rmse(true_reg, regs)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", default="data/processed/iedb_class1_human_nM.tsv")
    ap.add_argument("--out", default="results/phla_benchmark.json")
    ap.add_argument("--epochs-cnn", type=int, default=8)
    ap.add_argument("--epochs-gnn", type=int, default=8)
    ap.add_argument("--batch", type=int, default=512)
    ap.add_argument("--lr", type=float, default=1e-3)
    args = ap.parse_args()

    df = load_filtered_tsv(args.tsv)
    examples = aggregate(df)
    counts: dict[str, int] = defaultdict(int)
    for e in examples:
        counts[e.allele] += 1
    keep_alleles = {a for a, c in counts.items() if c >= MIN_PER_ALLELE}
    examples = [e for e in examples if e.allele in keep_alleles]
    alleles = sorted(keep_alleles)
    allele_map = {a: i for i, a in enumerate(alleles)}
    print(f"examples={len(examples)} alleles={len(alleles)}", flush=True)

    train, val, test = split_by_peptide(examples)
    print(f"train={len(train)} val={len(val)} test={len(test)}", flush=True)

    results: dict = {"n_examples": len(examples), "n_alleles": len(alleles),
                     "alleles": alleles,
                     "split": {"train": len(train), "val": len(val), "test": len(test)},
                     "binder_threshold_nm": BINDER_THRESHOLD_NM}

    # --- PSSM baseline (per allele) ---
    by_allele_train: dict[str, list[PHLAExample]] = defaultdict(list)
    for e in train:
        by_allele_train[e.allele].append(e)
    pssms = {a: AllelePSSM().fit([e.sequence for e in exs],
                                 np.array([e.log_ic50 for e in exs]))
             for a, exs in by_allele_train.items()}
    test_labels = np.array([e.binder for e in test])
    test_reg = np.array([e.log_ic50 for e in test])
    pssm_pred = np.array([-pssms[e.allele].predict([e.sequence])[0] for e in test])
    pssm_reg = np.array([pssms[e.allele].predict([e.sequence])[0] for e in test])
    results["pssm"] = {"auc": M.auc(test_labels, pssm_pred),
                       "auc0.1": M.auc_top_frac(test_labels, pssm_pred),
                       "ppv": M.ppv(test_labels, pssm_pred),
                       "srcc": M.srcc(test_reg, pssm_reg),
                       "rmse": M.rmse(test_reg, pssm_reg)}
    print("PSSM:", {k: round(v, 4) for k, v in results["pssm"].items()}, flush=True)

    # --- CNN ---
    cnn = PHLACNN(n_alleles=len(alleles))
    train_torch_model(cnn, train, val, allele_map, args.epochs_cnn, args.batch,
                      args.lr, with_graph=False)
    results["cnn"] = evaluate_torch(cnn, test, allele_map, with_graph=False)
    print("CNN:", {k: round(v, 4) for k, v in results["cnn"].items()}, flush=True)
    torch.save({"state": cnn.state_dict(), "alleles": alleles}, "results/phla_cnn.pt")

    # --- GNN ---
    gnn = PHLAGNN(n_alleles=len(alleles))
    train_torch_model(gnn, train, val, allele_map, args.epochs_gnn, args.batch,
                      args.lr, with_graph=True)
    results["gnn"] = evaluate_torch(gnn, test, allele_map, with_graph=True)
    print("GNN:", {k: round(v, 4) for k, v in results["gnn"].items()}, flush=True)
    torch.save({"state": gnn.state_dict(), "alleles": alleles}, "results/phla_gnn.pt")

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(results, indent=2))
    print(f"wrote {args.out}", flush=True)


if __name__ == "__main__":
    main()
