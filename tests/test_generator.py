"""Hermetic tests for the CPP generator and screening cascade."""
import torch

from peptidehlacpp.design.generator import (CPPGenerator, MAX_GEN_LEN, TOK,
                                            VOCAB, encode_seq,
                                            screen_candidates)


def test_encode_seq_bookends():
    toks = encode_seq("ACD")
    assert toks[0] == TOK["<s>"] and toks[-1] == TOK["<e>"]
    assert len(toks) == 5


def test_generator_forward_shape():
    m = CPPGenerator()
    logits = m(torch.tensor([encode_seq("ACDEFGHIK")[:-1]]))
    assert logits.shape == (1, len(encode_seq("ACDEFGHIK")) - 1, len(VOCAB))
    assert torch.isfinite(logits).all()


def test_screen_filters_charge_and_duplicates():
    cands = ["KKKKKKKKKK", "DDDDDDDDDD", "KKKKKKKKKK", "ACDEFGHIK"]
    passed = screen_candidates(
        cands, classifier_score=lambda s: 0.99,
        train_seqs={"ACDEFGHIK"}, min_score=0.5,
        charge_range=(2.0, 12.0), min_moment=0.0, max_hydropathy=10.0)
    seqs = [p["sequence"] for p in passed]
    assert seqs == ["KKKKKKKKKK"]  # acidic fails charge; training seq excluded; dup collapsed


def test_screen_respects_score_threshold():
    passed = screen_candidates(["KKKKKKKKK"], classifier_score=lambda s: 0.3,
                               train_seqs=set(), min_score=0.5)
    assert passed == []
