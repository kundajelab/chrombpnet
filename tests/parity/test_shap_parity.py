"""DeepSHAP parity against chrombpnet 1.x (kundajelab-shap TFDeepExplainer) goldens.

Goldens come from tests/goldens/make_goldens.py (legacy TF 2.12 container, CPU, content-seeded references):
    shap_<name>.npz           x (N, L, 4) int8, refs (N, K, L, 4) int8
    shap_<name>_<head>.npz    mult (N, K, L, 4) float32 (input-half multipliers), hyp (N, L, 4) float64
    fvals_<name>.npz          logits_x, logcounts_x, logits_refs, logcounts_refs

Environment:
    CHROMBPNET_GOLDENS         directory with the golden files (tests skip without it)
    CHROMBPNET_LEGACY_BIAS     the bias model <b> of the bias goldens (shap_bias_<b>*.npz, fvals_bias_<b>.npz; the
                               BIAS of run_goldens.sh); the bias tests skip without it
    CHROMBPNET_LEGACY_RESULTS  chrombpnet 1.x results tree with the real models and CHROMBPNET_LEGACY_PREFIX its run
                               name (see tests/conftest.py; model tests skip without them), or
    CHROMBPNET_MODEL_BIAS / CHROMBPNET_MODEL_NOBIAS   explicit model paths
    CHROMBPNET_PARITY_MAX_SEQS cap on sequences per model test (CPU runs of the 512x8 model are slow)
    CHROMBPNET_PARITY_BATCH_SEQS sequences per device step (default 8)
On GPUs also export XLA_FLAGS=--xla_gpu_deterministic_ops=true before running.
"""
import os

import numpy as np
import pytest

BIAS = os.environ.get("CHROMBPNET_LEGACY_BIAS")
BIAS_NAME = "bias_" + BIAS if BIAS else pytest.param("bias", marks=pytest.mark.skip(
    reason="CHROMBPNET_LEGACY_BIAS is not set"))
NAMES = (BIAS_NAME, "chrombpnet_nobias")
HEADS = ("counts", "profile")


def model_env(name):
    return "CHROMBPNET_MODEL_NOBIAS" if name == "chrombpnet_nobias" else "CHROMBPNET_MODEL_BIAS"


def model_rel(name, prefix):
    """Model path inside CHROMBPNET_LEGACY_RESULTS (the layout tests/goldens/run_goldens.sh reads)."""
    if name == "chrombpnet_nobias":
        return os.path.join("full_models", prefix, "models", "chrombpnet_nobias.h5")
    return os.path.join("bias_models", "bias_model_" + name[len("bias_"):], prefix, "models", prefix + "_bias.h5")


def golden(filename):
    root = os.environ.get("CHROMBPNET_GOLDENS")
    if not root:
        pytest.skip("CHROMBPNET_GOLDENS not set")
    path = os.path.join(root, filename)
    if not os.path.exists(path):
        pytest.skip("golden {} missing".format(path))
    return np.load(path)


def model_path(name):
    path = os.environ.get(model_env(name))
    if not path:
        root = os.environ.get("CHROMBPNET_LEGACY_RESULTS")
        prefix = os.environ.get("CHROMBPNET_LEGACY_PREFIX")
        if not root or not prefix:
            pytest.skip("CHROMBPNET_LEGACY_RESULTS + CHROMBPNET_LEGACY_PREFIX / {} not set".format(model_env(name)))
        path = os.path.join(root, model_rel(name, prefix))
    if not os.path.exists(path):
        pytest.skip("model {} missing".format(path))
    return path


def max_seqs(n):
    cap = os.environ.get("CHROMBPNET_PARITY_MAX_SEQS")
    return min(n, int(cap)) if cap else n


def rel_l2(a, b):
    a = np.asarray(a, np.float64)
    b = np.asarray(b, np.float64)
    return float(np.linalg.norm(a - b) / np.linalg.norm(b))


def softmax(a):
    e = np.exp(a - a.max(-1, keepdims=True))
    return e / e.sum(-1, keepdims=True)


def backend_tolerance():
    import jax
    return 1e-4 if jax.default_backend() == "cpu" else 1e-3


@pytest.mark.parametrize("name", NAMES)
def test_references_reproduce_goldens(name):
    from chrombpnet.evaluation.interpret.dinuc_shuffle import make_references
    d = golden("shap_{}.npz".format(name))
    refs = make_references(d["x"], num_shuffles=d["refs"].shape[1], seed=None)
    assert refs.dtype == d["refs"].dtype
    np.testing.assert_array_equal(refs, d["refs"])


@pytest.mark.parametrize("name", NAMES)
def test_t2_legacy_profile_multipliers_follow_midpoint_rule(name):
    """sum m*(x-r) of the LEGACY profile multipliers vs three readings of deep_tf.py."""
    d = golden("shap_{}.npz".format(name))
    g = golden("shap_{}_profile.npz".format(name))
    fv = golden("fvals_{}.npz".format(name))
    x, refs, mult = d["x"].astype(np.float64), d["refs"].astype(np.float64), g["mult"].astype(np.float64)
    lhs = np.sum(mult * (x[:, None] - refs), axis=(2, 3))
    lx = fv["logits_x"].astype(np.float64)
    lr = fv["logits_refs"].astype(np.float64)
    mnx = (lx - lx.mean(-1, keepdims=True))[:, None]
    mnr = lr - lr.mean(-1, keepdims=True)
    sx, sr = softmax(mnx), softmax(mnr)
    delta = mnx - mnr
    h1 = np.sum(np.where(np.abs(delta) < 1e-7, 0.0, 0.5 * (sx + sr)) * delta, axis=-1)  # midpoint
    h2 = np.sum(sx * delta, axis=-1)  # stop_gradient(softmax(x)) autodiff
    h3 = np.sum(sx * mnx, axis=-1) - np.sum(sr * mnr, axis=-1)  # full rescale T(x) - T(r)
    scale = np.maximum(np.abs(h1), 1.0)
    err = {k: float(np.median(np.abs(lhs - h) / scale)) for k, h in (("midpoint", h1), ("softmax_x", h2),
                                                                      ("full_rescale", h3))}
    print(name, err)
    assert err["midpoint"] < 1e-3, err
    if err["softmax_x"] > 1e-2 and err["full_rescale"] > 1e-2:
        assert err["midpoint"] * 10 < err["softmax_x"] and err["midpoint"] * 10 < err["full_rescale"], err
    else:
        # a barely trained model has near-flat profiles, where the three readings coincide
        print(name, "hypotheses not separable on this model")


@pytest.mark.parametrize("name", NAMES)
def test_legacy_counts_multipliers_sum_to_delta(name):
    d = golden("shap_{}.npz".format(name))
    g = golden("shap_{}_counts.npz".format(name))
    fv = golden("fvals_{}.npz".format(name))
    x, refs = d["x"].astype(np.float64), d["refs"].astype(np.float64)
    lhs = np.sum(g["mult"].astype(np.float64) * (x[:, None] - refs), axis=(2, 3))
    rhs = fv["logcounts_x"].astype(np.float64)[:, None, 0] - fv["logcounts_refs"].astype(np.float64)[..., 0]
    assert float(np.max(np.abs(lhs - rhs) / np.maximum(np.abs(rhs), 1.0))) < 1e-3


@pytest.fixture(scope="module")
def models():
    cache = {}

    def get(name):
        if name not in cache:
            from chrombpnet.training.utils.model_io import load_model
            cache[name] = load_model(model_path(name))
        return cache[name]
    return get


@pytest.mark.slow
@pytest.mark.parametrize("name", NAMES)
def test_model_predictions_match_legacy(name, models):
    d = golden("shap_{}.npz".format(name))
    fv = golden("fvals_{}.npz".format(name))
    n = max_seqs(len(d["x"]))
    import jax
    with jax.default_matmul_precision("highest"):
        logits, logcounts = models(name).predict(d["x"][:n].astype(np.float32), batch_size=16, verbose=0)
    atol = 1e-4 if jax.default_backend() == "cpu" else 1e-3
    np.testing.assert_allclose(logits, fv["logits_x"][:n], atol=atol, rtol=0)
    np.testing.assert_allclose(logcounts, fv["logcounts_x"][:n], atol=atol, rtol=0)


@pytest.mark.slow
@pytest.mark.parametrize("head", HEADS)
@pytest.mark.parametrize("name", NAMES)
def test_multipliers_and_hypothetical_match_legacy(name, head, models):
    from chrombpnet.evaluation.interpret.explainer import DeepLiftShap
    d = golden("shap_{}.npz".format(name))
    g = golden("shap_{}_{}.npz".format(name, head))
    n = max_seqs(len(d["x"]))
    x, refs = d["x"][:n], d["refs"][:n]
    batch_seqs = int(os.environ.get("CHROMBPNET_PARITY_BATCH_SEQS", 8))  # results do not depend on it
    explainer = DeepLiftShap(models(name), heads=[head], precision="highest", batch_seqs=batch_seqs)
    res = explainer.run(x, references=refs, return_multipliers=True)
    tol = backend_tolerance()
    mult_err = rel_l2(res.mult[head], g["mult"][:n])
    hyp_err = rel_l2(res.hyp[head], g["hyp"][:n])
    per_region = [rel_l2(res.hyp[head][i], g["hyp"][i]) for i in range(n)]
    print(name, head, "mult rel L2 {:.2e}  hyp rel L2 {:.2e}  worst region {:.2e}  additivity {}".format(
        mult_err, hyp_err, max(per_region), res.additivity[head]))
    assert mult_err <= tol
    assert hyp_err <= tol
