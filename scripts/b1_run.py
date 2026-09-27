"""B1 unseen-allele benchmark runner (PREREG_B1 + locked protocol details,
docs/PREREG_B1_PROTOCOL_DETAILS_2026-09-27.md). Resume-safe per phase:
each phase writes its own JSON/PT/NPZ checkpoint under results/b1_partial/
and skips if complete. Phases:
  arm3      pooled PSSM fit + held-out scoring (fast)
  train-cnn CNN training on all 52-allele rows (per-epoch checkpoints)
  arm1      mean-embedding CNN inference on held-out alleles
  embed     ESM-2 frozen embeddings (chunk checkpoints)
  arm2      logistic head on ESM-2 embeddings + held-out scoring
  finalize  falsifier + 9-mer dup audit + bootstrap + results/b1_unseen_allele.json
"""
import argparse, json, sys, time
from collections import defaultdict
from pathlib import Path
sys.path.insert(0, "src")
import numpy as np

from peptidehlacpp.data.iedb_dataset import aggregate, load_filtered_tsv
from peptidehlacpp.eval import metrics as M
from peptidehlacpp.models.pssm import AllelePSSM

PART = Path("results/b1_partial")
SEED_VAL, SEED_BOOT, SEED_SHUFFLE = 11, 23, 29
FALSIFIER_ALLELE = "HLA-A*32:07"


def load_data():
    df = load_filtered_tsv("data/processed/iedb_class1_human_nM.tsv")
    examples = aggregate(df)
    train_alleles = {r["allele"] for r in json.load(open("results/per_allele_analysis.json"))["per_allele"]}
    held = json.load(open("results/b1_allele_set.json"))["alleles"]
    held_alleles = [h["allele"] for h in held]
    train_ex = [e for e in examples if e.allele in train_alleles]
    by_allele = {a: [e for e in examples if e.allele == a] for a in held_alleles}
    return train_ex, by_allele, sorted(train_alleles)


def per_allele_metrics(labels, scores):
    labels = np.asarray(labels); scores = np.asarray(scores)
    n_pos = int(labels.sum()); n = len(labels)
    out = {"n": n, "n_pos": n_pos}
    if n_pos == 0 or n_pos == n:
        out.update(auroc=None, pauc01=None, calib_slope=None,
                   note="single-class allele: AUROC undefined (locked handling)")
        return out
    out["auroc"] = M.auc(labels, scores)
    out["pauc01"] = M.auc_top_frac(labels, scores, frac=0.1) if hasattr(M, "auc_top_frac") else None
    if n_pos >= 5 and n - n_pos >= 5:
        from sklearn.linear_model import LogisticRegression
        lr = LogisticRegression(max_iter=1000).fit(scores.reshape(-1, 1), labels)
        out["calib_slope"] = float(lr.coef_[0][0])
    else:
        out["calib_slope"] = None
    return out


def phase_arm3():
    out_p = PART / "arm3_per_allele.json"
    if out_p.exists():
        print("arm3 done, skip"); return
    train_ex, by_allele, _ = load_data()
    print(f"arm3: pooled PSSM on {len(train_ex)} rows", flush=True)
    pssm = AllelePSSM(ridge=1.0).fit([e.sequence for e in train_ex],
                                     np.array([e.log_ic50 for e in train_ex]))
    res = {}
    for a, exs in by_allele.items():
        scores = -pssm.predict([e.sequence for e in exs])
        m = per_allele_metrics([e.binder for e in exs], scores)
        m["scores"] = [round(float(s), 6) for s in scores]
        m["labels"] = [int(e.binder) for e in exs]
        res[a] = m
    json.dump(res, open(out_p, "w"))
    au = [v["auroc"] for v in res.values() if v["auroc"] is not None]
    print(f"arm3 mean AUROC {np.mean(au):.4f} over {len(au)} alleles", flush=True)


def phase_train_cnn():
    import torch
    from peptidehlacpp.models.cnn import PHLACNN
    from peptidehlacpp.training.train_phla import enc_batch
    train_ex, _, train_alleles = load_data()
    amap = {a: i for i, a in enumerate(train_alleles)}
    # 5% peptide-level val split from training rows, seed 11 (locked)
    rng = np.random.default_rng(SEED_VAL)
    peps = np.array(sorted({e.sequence for e in train_ex}))
    rng.shuffle(peps)
    val_peps = set(peps[: int(round(0.05 * len(peps)))])
    tr = [e for e in train_ex if e.sequence not in val_peps]
    va = [e for e in train_ex if e.sequence in val_peps]
    log_p = PART / "b1_cnn_train_log.json"
    log = json.load(open(log_p)) if log_p.exists() else {"epochs": [], "val_peps_seed": SEED_VAL}
    done = len(log["epochs"])
    if log.get("complete"):
        print("train-cnn done, skip"); return
    import torch.nn as nn
    torch.manual_seed(3)
    model = PHLACNN(n_alleles=len(train_alleles))
    if done:
        model.load_state_dict(torch.load(PART / f"b1_cnn_epoch{done}.pt")["state"])
    opt = torch.optim.AdamW(model.parameters(), lr=1e-3, weight_decay=1e-4)
    bce, mse = nn.BCEWithLogitsLoss(), nn.MSELoss()
    y_reg = torch.from_numpy(np.array([e.log_ic50 for e in tr], dtype=np.float32))
    y_cls = torch.from_numpy(np.array([e.binder for e in tr], dtype=np.float32))
    rng_t = np.random.default_rng(3 + done)
    best_auc, best_ep = (log.get("best_auc", -1), log.get("best_epoch", 0))
    n = len(tr)
    for ep in range(done, 8):
        model.train()
        order = rng_t.permutation(n)
        t0 = time.time()
        for s in range(0, n, 512):
            idx = order[s:s + 512]
            seqs = [tr[i].sequence for i in idx]
            als = np.array([tr[i].allele for i in idx])
            X, m, al = enc_batch(seqs, als, amap, False)
            pr, pc = model(X, m, al)
            loss = mse(pr, y_reg[idx]) + bce(pc, y_cls[idx])
            opt.zero_grad(); loss.backward(); opt.step()
        # val AUC on classifier logit
        model.eval(); scs = []
        with torch.no_grad():
            for s in range(0, len(va), 2048):
                chunk = va[s:s + 2048]
                X, m, al = enc_batch([e.sequence for e in chunk],
                                     np.array([e.allele for e in chunk]), amap, False)
                scs.append(model(X, m, al)[1].numpy())
        va_auc = M.auc(np.array([e.binder for e in va]), np.concatenate(scs))
        log["epochs"].append({"epoch": ep + 1, "val_auc": va_auc, "secs": round(time.time() - t0, 1)})
        torch.save({"state": model.state_dict()}, PART / f"b1_cnn_epoch{ep + 1}.pt")
        if va_auc > best_auc:
            best_auc, best_ep = va_auc, ep + 1
            torch.save({"state": model.state_dict()}, PART / "b1_cnn_best.pt")
        log["best_auc"], log["best_epoch"] = best_auc, best_ep
        json.dump(log, open(log_p, "w"), indent=1)
        print(f"epoch {ep+1}/8 val_auc={va_auc:.4f} ({time.time()-t0:.0f}s)", flush=True)
    log["complete"] = True
    json.dump(log, open(log_p, "w"), indent=1)


def phase_arm1():
    out_p = PART / "arm1_per_allele.json"
    if out_p.exists():
        print("arm1 done, skip"); return
    import torch
    from peptidehlacpp.models.cnn import PHLACNN
    from peptidehlacpp.training.train_phla import enc_batch
    _, by_allele, train_alleles = load_data()
    amap = {a: i for i, a in enumerate(train_alleles)}
    model = PHLACNN(n_alleles=len(train_alleles))
    model.load_state_dict(torch.load(PART / "b1_cnn_best.pt")["state"])
    model.eval()
    # mean of the 52 learned allele embeddings (locked zero-shot mean-embedding)
    with torch.no_grad():
        mean_vec = model.allele_emb.weight.mean(dim=0, keepdim=True)
    res = {}
    for a, exs in by_allele.items():
        scs = []
        with torch.no_grad():
            for s in range(0, len(exs), 2048):
                chunk = exs[s:s + 2048]
                from peptidehlacpp.features import stacked_enc, MAX_LEN_PHLA
                X = np.stack([stacked_enc(e.sequence, MAX_LEN_PHLA) for e in chunk])
                mask = np.zeros((len(chunk), MAX_LEN_PHLA), dtype=np.float32)
                for r, e in enumerate(chunk):
                    mask[r, : min(len(e.sequence), MAX_LEN_PHLA)] = 1.0
                Xt = torch.from_numpy(X); mt = torch.from_numpy(mask)
                # bypass allele embedding index: temporarily use index 0 then
                # substitute the pooled mean embedding (identical for all rows)
                alt = torch.zeros(len(chunk), dtype=torch.int64)
                h = Xt.transpose(1, 2)
                for b in model.blocks:
                    h = b(h)
                m1 = mt.unsqueeze(1)
                h = h * m1
                mean_pool = h.sum(dim=2) / m1.sum(dim=2).clamp(min=1.0)
                max_pool = h.masked_fill(m1 == 0, -1e4).amax(dim=2)
                pooled = torch.cat([mean_pool, max_pool,
                                    mean_vec.expand(len(chunk), -1)], dim=1)
                z = model.head(pooled)
                scs.append(model.cls_out(z).squeeze(-1).numpy())
        scores = np.concatenate(scs)
        m = per_allele_metrics([e.binder for e in exs], scores)
        m["scores"] = [round(float(s), 6) for s in scores]
        m["labels"] = [int(e.binder) for e in exs]
        res[a] = m
        print(f"arm1 {a} auroc={m['auroc']}", flush=True)
    json.dump(res, open(out_p, "w"))
    au = [v["auroc"] for v in res.values() if v["auroc"] is not None]
    print(f"arm1 mean AUROC {np.mean(au):.4f} over {len(au)} alleles", flush=True)


def phase_embed():
    """ESM-2 frozen mean-pooled embeddings for all unique peptides.
    Chunk-checkpointed: results/b1_partial/esm2_emb_partial.npz."""
    done_p = PART / "esm2_emb.npz"
    if done_p.exists():
        print("embed done, skip"); return
    import torch, esm
    train_ex, by_allele, _ = load_data()
    peps = sorted({e.sequence for e in train_ex} |
                  {e.sequence for exs in by_allele.values() for e in exs})
    print(f"embed: {len(peps)} unique peptides", flush=True)
    model, alphabet = esm.pretrained.esm2_t12_35M_UR50D()
    model.eval()
    bc = alphabet.get_batch_converter()
    part_p = PART / "esm2_emb_partial.npz"
    start = 0
    emb_chunks = []
    if part_p.exists():
        d = np.load(part_p)
        emb_chunks = [d["emb"]]
        start = int(d["done"])
        print(f"resume at {start}", flush=True)
    B = 64
    with torch.no_grad():
        for s in range(start, len(peps), B):
            chunk = peps[s:s + B]
            _, _, toks = bc([(f"p{i}", p) for i, p in enumerate(chunk)])
            out = model(toks, repr_layers=[12])
            rep = out["representations"][12]
            for r, p in enumerate(chunk):
                L = len(p)
                emb_chunks.append(rep[r, 1:L + 1].mean(0).numpy()[None, :])
            done_now = s + len(chunk)
            if (done_now // B) % 20 == 0 or done_now == len(peps):
                np.savez_compressed(part_p, emb=np.concatenate(emb_chunks), done=done_now)
                print(f"embed {done_now}/{len(peps)}", flush=True)
                emb_chunks = [np.concatenate(emb_chunks)]
    emb = np.concatenate(emb_chunks)
    np.savez_compressed(done_p, peps=np.array(peps), emb=emb)
    part_p.unlink(missing_ok=True)


def phase_arm2():
    out_p = PART / "arm2_per_allele.json"
    if out_p.exists():
        print("arm2 done, skip"); return
    from sklearn.linear_model import LogisticRegression
    train_ex, by_allele, _ = load_data()
    d = np.load(PART / "esm2_emb.npz", allow_pickle=True)
    peps = list(d["peps"]); emb = d["emb"]
    idx = {p: i for i, p in enumerate(peps)}
    Xtr = emb[[idx[e.sequence] for e in train_ex]]
    ytr = np.array([e.binder for e in train_ex], dtype=int)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(Xtr, ytr)
    res = {}
    for a, exs in by_allele.items():
        X = emb[[idx[e.sequence] for e in exs]]
        scores = clf.decision_function(X)
        m = per_allele_metrics([e.binder for e in exs], scores)
        m["scores"] = [round(float(s), 6) for s in scores]
        m["labels"] = [int(e.binder) for e in exs]
        res[a] = m
    json.dump(res, open(out_p, "w"))
    au = [v["auroc"] for v in res.values() if v["auroc"] is not None]
    print(f"arm2 mean AUROC {np.mean(au):.4f} over {len(au)} alleles", flush=True)


def phase_finalize():
    train_ex, by_allele, _ = load_data()
    arms = {k: json.load(open(PART / f"{k}_per_allele.json")) for k in ("arm1", "arm2", "arm3")}
    # falsifier: arm1 scores on FALSIFIER_ALLELE vs seed-29 label shuffle
    f = arms["arm1"][FALSIFIER_ALLELE]
    rng = np.random.default_rng(SEED_SHUFFLE)
    shuf = rng.permutation(np.array(f["labels"]))
    fals_auroc = M.auc(shuf, np.array(f["scores"]))
    falsifier_ok = 0.40 <= fals_auroc <= 0.60
    # 9-mer duplicate audit
    train9, train8 = set(), set()
    for e in train_ex:
        s = e.sequence
        if len(s) >= 9:
            for i in range(len(s) - 8):
                train9.add(s[i:i + 9])
        if len(s) == 8:
            train8.add(s)
    audit = {}
    for a, exs in by_allele.items():
        hits = 0
        for e in exs:
            s = e.sequence
            if len(s) >= 9:
                hit = any(s[i:i + 9] in train9 for i in range(len(s) - 8))
            else:
                hit = s in train8
            hits += hit
        audit[a] = {"overlap_frac": round(hits / len(exs), 4),
                    "flagged_gt_90pct": bool(hits / len(exs) > 0.90)}
    # primary endpoint + bootstrap over alleles
    rngb = np.random.default_rng(SEED_BOOT)
    summary = {}
    for arm, res in arms.items():
        aus = np.array([v["auroc"] for v in res.values() if v["auroc"] is not None])
        boots = [float(rngb.choice(aus, size=len(aus), replace=True).mean()) for _ in range(10000)]
        summary[arm] = {
            "mean_per_allele_auroc": float(aus.mean()),
            "ci95": [float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))],
            "n_alleles_scored": len(aus),
            "n_alleles_single_class": sum(1 for v in res.values() if v["auroc"] is None),
            "per_allele": {a: {k: v[k] for k in ("n", "n_pos", "auroc", "pauc01", "calib_slope")}
                           for a, v in res.items()},
        }
    def decision(m, lo):
        if m >= 0.80 and lo > 0.70:
            return "SUCCESS"
        if m >= 0.70:
            return "PARTIAL"
        return "FAILURE"
    out = {
        "prereg": ["docs/PREREG_B1_UNSEEN_ALLELE_2026-09-27.md",
                   "docs/PREREG_B1_PROTOCOL_DETAILS_2026-09-27.md"],
        "allele_set": "results/b1_allele_set.json (frozen 9e8a145)",
        "falsifier": {"allele": FALSIFIER_ALLELE, "seed": SEED_SHUFFLE,
                      "shuffled_auroc": fals_auroc, "in_0.40_0.60": falsifier_ok},
        "dup_audit_9mer": audit,
        "arms": summary,
        "decision": {arm: decision(s["mean_per_allele_auroc"], s["ci95"][0])
                     for arm, s in summary.items()},
        "g3_note": "G3 locked negative (0.6687 vs 0.9285 in-distribution) stands regardless of B1 outcome.",
    }
    json.dump(out, open("results/b1_unseen_allele.json", "w"), indent=1)
    print(json.dumps({k: out[k] for k in ("falsifier", "decision")}, indent=1))
    for arm, s in summary.items():
        print(arm, s["mean_per_allele_auroc"], s["ci95"])


def main():
    PART.mkdir(exist_ok=True)
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase", required=True,
                    choices=["arm3", "train-cnn", "arm1", "embed", "arm2", "finalize"])
    a = ap.parse_args()
    {"arm3": phase_arm3, "train-cnn": phase_train_cnn, "arm1": phase_arm1,
     "embed": phase_embed, "arm2": phase_arm2, "finalize": phase_finalize}[a.phase]()


if __name__ == "__main__":
    main()
