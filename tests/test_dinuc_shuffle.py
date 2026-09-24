"""Dinucleotide shuffle: bit-identical to deeplift 0.6.13, content-seeded references, shuffle properties."""
import zlib
from collections import Counter

import numpy as np
import pytest

from chrombpnet.evaluation.interpret import dinuc_shuffle as ds


# ---- verbatim copy of deeplift 0.6.13 deeplift/dinuc_shuffle.py (one-hot path), MIT License, Kundaje Lab ----
def _deeplift_one_hot_to_tokens(one_hot):
    tokens = np.tile(one_hot.shape[1], one_hot.shape[0])  # Vector of all D
    seq_inds, dim_inds = np.where(one_hot)
    tokens[seq_inds] = dim_inds
    return tokens


def _deeplift_tokens_to_one_hot(tokens, one_hot_dim):
    identity = np.identity(one_hot_dim + 1)[:, :-1]  # Last row is all 0s
    return identity[tokens]


def deeplift_dinuc_shuffle(seq, num_shufs=None, rng=None):
    seq_len, one_hot_dim = seq.shape
    arr = _deeplift_one_hot_to_tokens(seq)
    if not rng:
        rng = np.random.RandomState()
    chars, tokens = np.unique(arr, return_inverse=True)
    shuf_next_inds = []
    for t in range(len(chars)):
        mask = tokens[:-1] == t  # Excluding last char
        inds = np.where(mask)[0]
        shuf_next_inds.append(inds + 1)  # Add 1 for next token
    all_results = np.empty((num_shufs if num_shufs else 1, seq_len, one_hot_dim), dtype=seq.dtype)
    for i in range(num_shufs if num_shufs else 1):
        for t in range(len(chars)):
            inds = np.arange(len(shuf_next_inds[t]))
            inds[:-1] = rng.permutation(len(inds) - 1)  # Keep last index same
            shuf_next_inds[t] = shuf_next_inds[t][inds]
        counters = [0] * len(chars)
        ind = 0
        result = np.empty_like(tokens)
        result[0] = tokens[ind]
        for j in range(1, len(tokens)):
            t = tokens[ind]
            ind = shuf_next_inds[t][counters[t]]
            counters[t] += 1
            result[j] = tokens[ind]
        all_results[i] = _deeplift_tokens_to_one_hot(chars[result], one_hot_dim)
    return all_results if num_shufs else all_results[0]
# ---- end of copy ----


def content_seed(onehot):
    """tests/goldens/make_goldens.py content_seed."""
    return zlib.crc32(np.ascontiguousarray(onehot, dtype=np.int8).tobytes()) & 0xFFFFFFFF


def golden_make_refs(seqs, num_shufs=20):
    """tests/goldens/make_goldens.py make_refs, with the deeplift copy above."""
    refs = np.empty((seqs.shape[0], num_shufs) + seqs.shape[1:], dtype=np.int8)
    for i, s in enumerate(seqs):
        refs[i] = deeplift_dinuc_shuffle(s, num_shufs=num_shufs, rng=np.random.RandomState(content_seed(s)))
    return refs


def random_onehot(rng, n, length, n_frac=0.0):
    tokens = rng.randint(0, 4, (n, length))
    if n_frac:
        tokens[rng.rand(n, length) < n_frac] = 4
    return np.eye(5, dtype=np.int8)[:, :4][tokens]


def with_n_runs(rng, length):
    x = random_onehot(rng, 1, length)[0]
    x[10:40] = 0
    x[length // 2:length // 2 + 5] = 0
    x[-3:] = 0
    return x


def tokens_of(onehot):
    return ds.one_hot_to_tokens(onehot)


def dinucs(tokens):
    return Counter(zip(tokens[:-1].tolist(), tokens[1:].tolist()))


@pytest.fixture(scope="module")
def seqs():
    rng = np.random.RandomState(0)
    x = random_onehot(rng, 6, 2114, n_frac=0.02)
    x[0] = with_n_runs(rng, 2114)
    x[1, :500] = np.eye(4, dtype=np.int8)[1]  # homopolymer run
    return x


@pytest.mark.parametrize("use_numba", [True, False])
def test_bit_identical_to_deeplift(seqs, use_numba):
    for s in seqs:
        for seed in (0, 1, 7, 1234, 2**31 - 1):
            ours = ds.dinuc_shuffle(s, num_shufs=20, rng=np.random.RandomState(seed), use_numba=use_numba)
            ref = deeplift_dinuc_shuffle(s, num_shufs=20, rng=np.random.RandomState(seed))
            assert ours.dtype == np.int8
            np.testing.assert_array_equal(ours, ref)


def test_single_shuffle_and_string_paths():
    rng = np.random.RandomState(3)
    s = random_onehot(rng, 1, 300)[0]
    one = ds.dinuc_shuffle(s, rng=np.random.RandomState(5))
    assert one.shape == s.shape
    np.testing.assert_array_equal(one, deeplift_dinuc_shuffle(s, rng=np.random.RandomState(5)))
    shufs = ds.dinuc_shuffle("ACGTTGCAAC" * 20, num_shufs=3, rng=np.random.RandomState(1))
    assert len(shufs) == 3 and all(len(t) == 200 and t[0] == "A" for t in shufs)


def test_numba_walk_equals_python_walk(seqs):
    for s in seqs:
        a = ds.dinuc_shuffle(s, num_shufs=5, rng=np.random.RandomState(11), use_numba=True)
        b = ds.dinuc_shuffle(s, num_shufs=5, rng=np.random.RandomState(11), use_numba=False)
        np.testing.assert_array_equal(a, b)


def test_make_references_matches_golden_convention(seqs):
    np.testing.assert_array_equal(ds.make_references(seqs, seed=None), golden_make_refs(seqs))


def test_make_references_seeding(seqs):
    refs = ds.make_references(seqs, num_shuffles=20, seed=1234)
    assert refs.shape == (len(seqs), 20) + seqs.shape[1:] and refs.dtype == np.int8
    for s, r in zip(seqs[:2], refs[:2]):
        rng = np.random.RandomState((content_seed(s) ^ 1234) & 0xFFFFFFFF)
        np.testing.assert_array_equal(r, deeplift_dinuc_shuffle(s, num_shufs=20, rng=rng))
    assert not np.array_equal(refs, ds.make_references(seqs, seed=4321))
    np.testing.assert_array_equal(refs, ds.make_references(seqs, seed=1234))


def test_make_references_independent_of_order_and_batching(seqs):
    refs = ds.make_references(seqs)
    order = np.array([3, 0, 5, 1, 4, 2])
    np.testing.assert_array_equal(ds.make_references(seqs[order]), refs[order])
    parts = np.concatenate([ds.make_references(seqs[:2]), ds.make_references(seqs[2:])])
    np.testing.assert_array_equal(parts, refs)


def test_make_references_casts_to_int8(seqs):
    refs = ds.make_references(seqs[:2].astype(np.float32), num_shuffles=3)
    assert refs.dtype == np.int8
    np.testing.assert_array_equal(refs, ds.make_references(seqs[:2], num_shuffles=3))


def test_shuffle_properties(seqs):
    refs = ds.make_references(seqs, num_shuffles=20)
    for s, rs in zip(seqs, refs):
        t = tokens_of(s)
        for r in rs:
            assert set(np.unique(r)) <= {0, 1}
            assert (r.sum(-1) <= 1).all()
            u = tokens_of(r)
            assert u[0] == t[0]
            assert dinucs(u) == dinucs(t)
            assert Counter(u.tolist()) == Counter(t.tolist())  # includes the number of N positions
        assert not all(np.array_equal(r, s) for r in rs)


def test_n_positions_are_shuffled_as_a_fifth_character():
    rng = np.random.RandomState(8)
    s = with_n_runs(rng, 400)
    for r in ds.dinuc_shuffle(s, num_shufs=10, rng=np.random.RandomState(2)):
        n_mask = r.sum(-1) == 0
        assert n_mask.sum() == (s.sum(-1) == 0).sum()
        assert dinucs(tokens_of(r)) == dinucs(tokens_of(s))


def test_periodic_sequence_shuffles_to_itself():
    s = np.eye(4, dtype=np.int8)[np.tile([0, 1], 1057)]  # (AC)^1057
    refs = ds.make_references(s[None], num_shuffles=20)
    for r in refs[0]:
        np.testing.assert_array_equal(r, s)
