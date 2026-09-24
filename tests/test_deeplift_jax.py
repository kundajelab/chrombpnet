"""JAX DeepLIFT/DeepSHAP rules on tiny BPNet-shaped Keras 3 models (CPU).

Checks the rescale rule, the legacy profile weighting, summation-to-delta, the paired vs full (TFDeepExplainer)
layout, an independent NumPy float64 brute-force DeepLIFT, and the hypothetical projection.
"""
import warnings

import chrombpnet  # noqa: F401  (KERAS_BACKEND=jax)
import jax
import jax.numpy as jnp
import keras
import numpy as np
import pytest
from keras import layers

from chrombpnet.evaluation.interpret import deeplift_jax as dl
from chrombpnet.evaluation.interpret import shap_utils
from chrombpnet.evaluation.interpret.dinuc_shuffle import make_references
from chrombpnet.evaluation.interpret.explainer import DeepLiftShap, auto_batch_seqs, check_model
from chrombpnet.training.utils.layers import LogSumExp

L, F, ND, OUTLEN, K = 256, 8, 2, 100, 4


def build_bpnet(inputlen=L, filters=F, n_dil=ND, outlen=OUTLEN, seed=0, prefix="", dtype=None, name=None):
    """Same layer structure as chrombpnet's bpnet_model, small."""
    keras.utils.set_random_seed(seed)
    kw = {"dtype": dtype} if dtype else {}
    inp = keras.Input((inputlen, 4), name="sequence")
    x = layers.Conv1D(filters, 21, activation="relu", name=prefix + "bpnet_1st_conv", **kw)(inp)
    for i in range(1, n_dil + 1):
        conv_x = layers.Conv1D(filters, 3, activation="relu", dilation_rate=2 ** i, name=prefix + "bpnet_{}conv".format(i),
                               **kw)(x)
        crop = (x.shape[1] - conv_x.shape[1]) // 2
        x = layers.Cropping1D(crop, name=prefix + "bpnet_{}crop".format(i), **kw)(x)
        x = layers.add([conv_x, x], **kw)
    prof = layers.Conv1D(1, 75, name=prefix + "prof_out_precrop", **kw)(x)
    cropsize = prof.shape[1] // 2 - outlen // 2
    prof = layers.Cropping1D(cropsize, name=prefix + "logits_profile_predictions_preflatten", **kw)(prof)
    profile_out = layers.Flatten(name=prefix + "logits_profile_predictions", **kw)(prof)
    gap = layers.GlobalAveragePooling1D(name=prefix + "gap", **kw)(x)
    count_out = layers.Dense(1, name=prefix + "logcount_predictions", **kw)(gap)
    model = keras.Model(inp, [profile_out, count_out], name=name)
    # non-zero biases so ReLUs switch inside the sequence
    rng = np.random.RandomState(seed + 100)
    model.set_weights([w if w.ndim > 1 else rng.normal(0, 0.05, w.shape).astype(w.dtype)
                       for w in model.get_weights()])
    return model


def onehot(tokens):
    return np.eye(4, dtype=np.int8)[tokens]


@pytest.fixture(scope="module")
def model():
    return build_bpnet()


@pytest.fixture(scope="module")
def data():
    rng = np.random.RandomState(1)
    x = onehot(rng.randint(0, 4, (5, L)))
    refs = make_references(x, num_shuffles=K, seed=1234)
    return x, refs


@pytest.fixture(scope="module")
def result(model, data):
    x, refs = data
    return DeepLiftShap(model, batch_seqs=2).run(x, references=refs, return_multipliers=True)


def predict(model, x):
    logits, logcounts = model(np.asarray(x, np.float32), training=False)
    return np.asarray(logits, np.float64), np.asarray(logcounts, np.float64)[:, 0]


def rel_l2(a, b):
    return float(np.linalg.norm(np.asarray(a, np.float64) - b) / np.linalg.norm(b))


# ------------------------------------------------------------------ NumPy float64 brute-force DeepLIFT
def np_conv(x, w, b, d):
    k = w.shape[0]
    lout = x.shape[1] - d * (k - 1)
    out = np.zeros((x.shape[0], lout, w.shape[2]))
    for j in range(k):
        out += x[:, j * d:j * d + lout] @ w[j]
    return out + b


def np_conv_t(g, w, d, lin):
    dx = np.zeros((g.shape[0], lin, w.shape[1]))
    for j in range(w.shape[0]):
        dx[:, j * d:j * d + g.shape[1]] += g @ w[j].T
    return dx


def np_weights(model, prefix=""):
    get = lambda n: [np.asarray(v, np.float64) for v in model.get_layer(prefix + n).get_weights()]
    convs = [get("bpnet_1st_conv") + [1]] + [get("bpnet_{}conv".format(i)) + [2 ** i] for i in range(1, ND + 1)]
    return {"convs": convs, "prof": get("prof_out_precrop"), "dense": get("logcount_predictions")}


def np_forward(wts, x):
    zs, lens = [], []
    (w1, b1, _) = wts["convs"][0]
    z = np_conv(x, w1, b1, 1)
    zs.append(z)
    h = np.maximum(z, 0)
    for (w, b, d) in wts["convs"][1:]:
        lens.append(h.shape[1])
        z = np_conv(h, w, b, d)
        zs.append(z)
        c = np.maximum(z, 0)
        crop = (h.shape[1] - c.shape[1]) // 2
        h = c + h[:, crop:h.shape[1] - crop]
    wp, bp = wts["prof"]
    p = np_conv(h, wp, bp, 1)
    cs = p.shape[1] // 2 - OUTLEN // 2
    logits = p[:, cs:cs + OUTLEN, 0]
    wd, bd = wts["dense"]
    logcounts = (h.mean(1) @ wd + bd)[:, 0]
    return logits, logcounts, {"zs": zs, "lens": lens, "h_len": h.shape[1], "p_len": p.shape[1], "cs": cs,
                               "in_len": x.shape[1]}


def np_rescale(zx, zr):
    d = zx - zr
    small = np.abs(d) < dl.EPS_RELU
    return np.where(small, (zx > 0).astype(float), (np.maximum(zx, 0) - np.maximum(zr, 0)) / np.where(small, 1, d))


def np_deeplift(wts, x, r, head):
    """Multipliers of the input rows of pairs (x[i], r[i]) by explicit backprop with the rescale rule."""
    lx, cx, cache_x = np_forward(wts, x)
    lr, cr, cache_r = np_forward(wts, r)
    ms = [np_rescale(zx, zr) for zx, zr in zip(cache_x["zs"], cache_r["zs"])]
    n = x.shape[0]
    if head == "counts":
        g_logits, g_counts = np.zeros_like(lx), np.ones((n, 1))
    else:
        mnx = lx - lx.mean(1, keepdims=True)
        mnr = lr - lr.mean(1, keepdims=True)
        sx = np.exp(mnx) / np.exp(mnx).sum(1, keepdims=True)
        sr = np.exp(mnr) / np.exp(mnr).sum(1, keepdims=True)
        w = np.where(np.abs(mnx - mnr) < dl.EPS_MUL, 0.0, 0.5 * (sx + sr))
        g_logits, g_counts = w - w.mean(1, keepdims=True), np.zeros((n, 1))
    wp, _ = wts["prof"]
    gp = np.zeros((n, cache_x["p_len"], 1))
    gp[:, cache_x["cs"]:cache_x["cs"] + OUTLEN, 0] = g_logits
    gh = np_conv_t(gp, wp, 1, cache_x["h_len"])
    wd, _ = wts["dense"]
    gh += (g_counts @ wd.T)[:, None, :] / cache_x["h_len"]
    for i in range(len(wts["convs"]) - 1, 0, -1):
        w, _, d = wts["convs"][i]
        lin = cache_x["lens"][i - 1]
        gz = gh * ms[i]
        g_in = np_conv_t(gz, w, d, lin)
        crop = (lin - gz.shape[1]) // 2
        g_in[:, crop:lin - crop] += gh
        gh = g_in
    w1, _, _ = wts["convs"][0]
    return np_conv_t(gh * ms[0], w1, 1, cache_x["in_len"])


# ------------------------------------------------------------------ rules
def test_rescale_rule_cases():
    zx = jnp.array([1.0, -1.0, 0.5, -0.5, 0.3, 2e-7, -4e-7, 3.0], jnp.float32)
    zr = jnp.array([2.0, -2.0, -0.5, 0.5, 0.3, -3e-7, 1e-7, 3.0 - 2e-6], jnp.float32)
    z = jnp.concatenate([zx, zr])
    g = jax.grad(lambda z: jnp.sum(dl.dl_relu(z)))(z)
    expect_x = np.array([1.0, 0.0, 0.5, 0.5, 1.0, 1.0, 0.0, 1.0], np.float32)
    np.testing.assert_allclose(np.asarray(g[:8]), expect_x, rtol=1e-6)
    np.testing.assert_array_equal(np.asarray(g[8:]), 0.0)
    g_full = jax.grad(lambda z: jnp.sum(dl.dl_relu_full(z)))(z)
    np.testing.assert_array_equal(np.asarray(g_full[:8]), np.asarray(g[:8]))
    # reference half: same secant, its own fallback (zr > 0)
    expect_r = np.array([1.0, 0.0, 0.5, 0.5, 1.0, 0.0, 1.0, 1.0], np.float32)
    np.testing.assert_allclose(np.asarray(g_full[8:]), expect_r, rtol=1e-6)
    np.testing.assert_array_equal(np.asarray(dl.dl_relu(z)), np.maximum(np.asarray(z), 0))


def test_profile_weights_midpoint_rule():
    rng = np.random.RandomState(0)
    mn_x = jnp.asarray(rng.normal(size=(3, 50)), jnp.float32)
    mn_r = mn_x.at[:, :5].add(jnp.asarray(rng.normal(size=(3, 5)), jnp.float32))
    mn_r = mn_r.at[:, 10].set(mn_x[:, 10] + 5e-8)  # |delta| < 1e-7 -> zero weight
    s_x, s_r = jax.nn.softmax(mn_x, -1), jax.nn.softmax(mn_r, -1)
    w = dl.profile_weights(mn_x, mn_r, "chrombpnet")
    expect = np.where(np.abs(np.asarray(mn_x - mn_r)) < 1e-7, 0.0, 0.5 * np.asarray(s_x + s_r))
    np.testing.assert_allclose(np.asarray(w), expect, rtol=1e-6)
    assert np.all(np.asarray(w)[:, 10] == 0)
    np.testing.assert_allclose(np.asarray(dl.profile_weights(mn_x, mn_r, "chrombpnet_tf")), expect,
                               rtol=1e-3, atol=1e-6)
    np.testing.assert_array_equal(np.asarray(dl.profile_weights(mn_x, mn_r, "softmax_x")), np.asarray(s_x))
    with pytest.raises(ValueError):
        dl.profile_weights(mn_x, mn_r, "nope")


# ------------------------------------------------------------------ explainer
def test_counts_summation_to_delta(model, data, result):
    x, refs = data
    _, cx = predict(model, x)
    _, cr = predict(model, refs.reshape(-1, L, 4))
    cr = cr.reshape(len(x), K)
    lhs = np.sum(result.mult["counts"] * (x[:, None] - refs), axis=(2, 3))
    rhs = cx[:, None] - cr
    np.testing.assert_allclose(lhs, rhs, rtol=1e-4, atol=1e-4 * np.abs(rhs).max())
    assert result.additivity["counts"]["max_rel_err"] < 1e-4


def test_profile_linearized_identity(model, data, result):
    x, refs = data
    lx, _ = predict(model, x)
    lr, _ = predict(model, refs.reshape(-1, L, 4))
    mnx = lx - lx.mean(1, keepdims=True)
    mnr = (lr - lr.mean(1, keepdims=True)).reshape(len(x), K, -1)
    w = np.asarray(dl.profile_weights(jnp.asarray(np.repeat(mnx, K, 0)), jnp.asarray(mnr.reshape(len(x) * K, -1))))
    rhs = np.sum(w.reshape(len(x), K, -1) * (mnx[:, None] - mnr), axis=-1)
    lhs = np.sum(result.mult["profile"] * (x[:, None] - refs), axis=(2, 3))
    np.testing.assert_allclose(lhs, rhs, rtol=1e-4, atol=1e-4 * np.abs(rhs).max())
    # it is the linearised delta, not T(x) - T(r): the dropped softmax share is visible
    t = lambda mn: np.sum(np.exp(mn) / np.exp(mn).sum(-1, keepdims=True) * mn, axis=-1)
    assert np.abs(lhs - (t(mnx)[:, None] - t(mnr))).max() > 100 * np.abs(lhs - rhs).max()


def test_paired_equals_full_variant(model, data, result):
    x, refs = data
    full = DeepLiftShap(model, batch_seqs=5, variant="full").run(x, references=refs, return_multipliers=True)
    for head in ("counts", "profile"):
        np.testing.assert_allclose(full.mult[head], result.mult[head], rtol=1e-5,
                                   atol=1e-6 * np.abs(result.mult[head]).max())


@pytest.mark.parametrize("head", ["counts", "profile"])
def test_matches_numpy_bruteforce_deeplift(model, data, result, head):
    x, refs = data
    wts = np_weights(model)
    xs = np.repeat(x, K, 0).astype(np.float64)
    rs = refs.reshape(-1, L, 4).astype(np.float64)
    oracle = np_deeplift(wts, xs, rs, head).reshape(len(x), K, L, 4)
    assert rel_l2(result.mult[head], oracle) < 1e-5
    lx, cx, _ = np_forward(wts, xs[:3])
    klx, kcx = predict(model, xs[:3])
    np.testing.assert_allclose(klx, lx, rtol=1e-5, atol=1e-6)
    np.testing.assert_allclose(kcx, cx, rtol=1e-5, atol=1e-6)


def test_rescale_rule_differs_from_plain_gradient(model, data, result):
    x, _ = data
    plain = jax.grad(lambda v: jnp.sum(model(v, training=False)[1]))(jnp.asarray(x, jnp.float32))
    diff = np.abs(result.mult["counts"] - np.asarray(plain)[:, None]).max()
    assert diff > 1e-3 * np.abs(np.asarray(plain)).max()


def test_hypothetical_matches_combine_mult_and_diffref(data, result):
    x, refs = data
    for head in ("counts", "profile"):
        oracle = np.stack([shap_utils.combine_mult_and_diffref([result.mult[head][i]], [x[i].astype(float)],
                                                               [refs[i]])[0] for i in range(len(x))])
        np.testing.assert_allclose(result.hyp[head], oracle, rtol=1e-5, atol=1e-7 * np.abs(oracle).max())


def test_periodic_sequence_has_zero_contributions(model):
    x = onehot(np.tile([0, 1], L // 2))[None]  # (AC)^n: every dinucleotide shuffle equals the input
    explainer = DeepLiftShap(model, batch_seqs=1, num_shuffles=K)
    res = explainer.run(x, return_multipliers=True, return_references=True)
    assert np.all(res.refs == x[:, None])
    for head in ("counts", "profile"):
        assert np.all(x * res.hyp[head] == 0)
    # every delta is 0: counts multipliers fall back to the plain gradient, profile weights are zeroed
    plain = jax.grad(lambda v: jnp.sum(model(v, training=False)[1]))(jnp.asarray(x, jnp.float32))
    np.testing.assert_allclose(res.mult["counts"], np.broadcast_to(np.asarray(plain)[:, None], res.mult["counts"].shape),
                               rtol=1e-5, atol=1e-7)
    assert np.all(res.mult["profile"] == 0)


def test_batching_and_padding_invariance(model, data, result):
    x, refs = data
    for bs in (1, 3, None):
        explainer = DeepLiftShap(model, batch_seqs=bs)
        hyp = explainer.run(x, references=refs).hyp
        assert len(explainer._steps) == 1
        for head in ("counts", "profile"):
            np.testing.assert_allclose(hyp[head], result.hyp[head], rtol=1e-5,
                                       atol=1e-6 * np.abs(result.hyp[head]).max())


def test_generated_references_are_content_seeded(model, data):
    x, _ = data
    explainer = DeepLiftShap(model, heads=["counts"], batch_seqs=2, num_shuffles=K)
    a = explainer.run(x, seed=7, return_references=True)
    np.testing.assert_array_equal(a.refs, make_references(x, num_shuffles=K, seed=7))
    b = explainer.run(x[::-1], seed=7)
    np.testing.assert_allclose(b.hyp["counts"][::-1], a.hyp["counts"], rtol=1e-6, atol=1e-8)


def test_activations_restored_and_model_unchanged(model, data):
    x, _ = data
    before = predict(model, x)
    DeepLiftShap(model, batch_seqs=5, num_shuffles=2).explain(x)
    for layer in model.layers:
        if isinstance(layer, layers.Conv1D):
            assert keras.activations.serialize(layer.activation) in ("relu", "linear")
    after = predict(model, x)
    np.testing.assert_array_equal(before[0], after[0])


def test_profile_weighting_alternatives(model, data):
    x, refs = data
    lx, _ = predict(model, x)
    lr, _ = predict(model, refs.reshape(-1, L, 4))
    mnx = np.repeat(lx - lx.mean(1, keepdims=True), K, 0)
    mnr = lr - lr.mean(1, keepdims=True)
    for mode in ("softmax_x", "tangermeme", "chrombpnet_tf"):
        mult = DeepLiftShap(model, heads=["profile"], profile_weighting=mode, batch_seqs=5).multipliers(
            x, references=refs)["profile"]
        lhs = np.sum(mult * (x[:, None] - refs), axis=(2, 3)).reshape(-1)
        w = np.asarray(dl.profile_weights(jnp.asarray(mnx, jnp.float32), jnp.asarray(mnr, jnp.float32), mode))
        rhs = np.sum(w * (mnx - mnr), axis=-1)
        np.testing.assert_allclose(lhs, rhs, rtol=1e-4, atol=1e-4 * np.abs(rhs).max())


def test_additivity_warning(model, data):
    x, refs = data
    explainer = DeepLiftShap(model, heads=["counts"], batch_seqs=5, additivity_threshold=1e-12)
    with pytest.warns(UserWarning, match="summation-to-delta"):
        explainer.run(x, references=refs)
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        DeepLiftShap(model, heads=["counts"], batch_seqs=5).run(x, references=refs)


def test_mixed_precision_model_is_explained_in_float32(model, data, result):
    x, refs = data
    bf16 = build_bpnet(dtype="mixed_bfloat16")
    bf16.set_weights(model.get_weights())
    assert bf16.get_layer("bpnet_1st_conv").compute_dtype == "bfloat16"
    res = DeepLiftShap(bf16, batch_seqs=2).run(x, references=refs)
    assert bf16.get_layer("bpnet_1st_conv").compute_dtype == "bfloat16"
    for head in ("counts", "profile"):
        np.testing.assert_allclose(res.hyp[head], result.hyp[head], rtol=1e-5,
                                   atol=1e-6 * np.abs(result.hyp[head]).max())


# ------------------------------------------------------------------ model support
def build_full_model():
    """chrombpnet_with_bias_model structure: no-bias + bias submodels, Add profile, logsumexp counts."""
    nobias = build_bpnet(seed=3, prefix="wo_bias_", name="model_wo_bias")
    bias = build_bpnet(seed=4, filters=4, name="model_bias")
    inp = keras.Input((L, 4), name="sequence")
    b_out, nb_out = bias(inp), nobias(inp)
    profile = layers.Add(name="logits_profile_predictions")([nb_out[0], b_out[0]])
    counts = LogSumExp(name="logcount_predictions")(layers.Concatenate(axis=-1)([nb_out[1], b_out[1]]))
    return keras.Model(inp, [profile, counts]), nobias, bias


def test_full_model_profile_supported_counts_refused(data):
    x, refs = data
    full, nobias, bias = build_full_model()
    with pytest.raises(NotImplementedError, match="_nobias.h5"):
        DeepLiftShap(full, heads=["counts", "profile"])
    res = DeepLiftShap(full, heads=["profile"], batch_seqs=5).run(x, references=refs, return_multipliers=True)
    assert res.additivity["profile"]["max_rel_err"] < 1e-4
    # the swap reaches the nested submodels: summation-to-delta of the full-model profile head
    lx = np.asarray(full(x.astype(np.float32))[0], np.float64)
    lr = np.asarray(full(refs.reshape(-1, L, 4).astype(np.float32))[0], np.float64)
    mnx = np.repeat(lx - lx.mean(1, keepdims=True), K, 0)
    mnr = lr - lr.mean(1, keepdims=True)
    w = np.asarray(dl.profile_weights(jnp.asarray(mnx), jnp.asarray(mnr)))
    rhs = np.sum(w * (mnx - mnr), axis=-1)
    lhs = np.sum(res.mult["profile"] * (x[:, None] - refs), axis=(2, 3)).reshape(-1)
    np.testing.assert_allclose(lhs, rhs, rtol=1e-4, atol=1e-4 * np.abs(rhs).max())
    DeepLiftShap(nobias, heads=["counts", "profile"])  # the no-bias submodel is fine


def test_unsupported_layers_raise():
    inp = keras.Input((64, 4))
    h = layers.Conv1D(4, 5, activation="sigmoid")(inp)
    out = layers.Dense(1)(layers.GlobalAveragePooling1D()(h))
    with pytest.raises(NotImplementedError, match="sigmoid"):
        check_model(keras.Model(inp, [layers.Flatten()(h), out]), ["counts"])
    h = layers.BatchNormalization()(layers.Conv1D(4, 5, activation="relu")(inp))
    out = layers.Dense(1)(layers.GlobalAveragePooling1D()(h))
    with pytest.raises(NotImplementedError, match="BatchNormalization"):
        check_model(keras.Model(inp, [layers.Flatten()(h), out]), ["profile"])


def test_auto_batch_seqs(model):
    # small models are capped at 32 sequences per step
    assert auto_batch_seqs(model) == 32
    # chrombpnet_nobias-sized (512 filters, 8 dilated layers, 2114 bp): ~4 GB per sequence (estimate) in 12 GB
    wide = build_bpnet(filters=512, n_dil=8, inputlen=2114)
    assert auto_batch_seqs(wide) == 2
    assert auto_batch_seqs(wide, n_heads=2) == 1

    class SmallDevice:
        def memory_stats(self):
            return {"bytes_limit": 3e9, "bytes_in_use": 1e9}

    # only 60% of the 2 GB still free may be used
    assert auto_batch_seqs(wide, device=SmallDevice()) == 1


def test_resolve_precision():
    from chrombpnet.evaluation.interpret.explainer import resolve_precision
    assert resolve_precision("auto", backend="cpu") == "highest"
    assert resolve_precision("auto", backend="gpu") == "default"
    assert resolve_precision("highest", backend="gpu") == "highest"