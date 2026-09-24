"""Keras 3 / JAX predictions of legacy (TF-Keras 2.x) chrombpnet 1.x models against goldens from the legacy stack.

Goldens: the directory written by `tests/goldens/make_goldens.py predict` (pred_inputs.npz + pred_<name>.npz, CPU,
float32), given by CHROMBPNET_GOLDENS. Models: the chrombpnet 1.x results tree given by CHROMBPNET_LEGACY_RESULTS,
with run name CHROMBPNET_LEGACY_PREFIX (see tests/conftest.py): pred_bias_<b>.npz is checked against
bias_models/bias_model_<b>/<prefix>/models/<prefix>_bias.h5 and pred_<name>.npz for the full models against
full_models/<prefix>/models/<name>.h5. Skipped unless all three are set. Run e.g.
`pytest -rP tests/parity/test_predict_parity.py` to see the observed maximum differences.
"""
import glob
import os

import numpy as np
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND before keras is imported)
import jax
from scipy.special import logsumexp

TOL = 1e-4  # max abs difference of logits / logcounts against CPU goldens at matmul precision 'highest'
FULL_MODELS = ["bias_model_scaled", "chrombpnet", "chrombpnet_nobias"]
LEGACY_BIAS = os.environ.get("CHROMBPNET_LEGACY_BIAS")
# parametrization when CHROMBPNET_GOLDENS is unset (every test then skips)
DEFAULT_NAMES = (["bias_" + LEGACY_BIAS] if LEGACY_BIAS else []) + FULL_MODELS


def env_dir(var):
    path = os.environ.get(var)
    if not path or not os.path.isdir(path):
        pytest.skip("{} is not set or is not a directory".format(var))
    return path


def golden_names():
    path = os.environ.get("CHROMBPNET_GOLDENS")
    if not path or not os.path.isdir(path):
        return DEFAULT_NAMES
    files = glob.glob(os.path.join(path, "pred_*.npz"))
    names = sorted(os.path.basename(f)[len("pred_"):-len(".npz")] for f in files)
    return [n for n in names if n != "inputs"] or DEFAULT_NAMES


def model_path(root, prefix, name):
    if name in FULL_MODELS:
        return os.path.join(root, "full_models", prefix, "models", name + ".h5")
    assert name.startswith("bias_"), name
    return os.path.join(root, "bias_models", "bias_model_" + name[len("bias_"):], prefix, "models",
                        prefix + "_bias.h5")


@pytest.fixture(scope="module")
def goldens():
    return env_dir("CHROMBPNET_GOLDENS")


@pytest.fixture(scope="module")
def predict(legacy_results, legacy_prefix):
    from chrombpnet.training.utils import model_io
    cache = {}

    def run(name, seqs):
        key = (name, len(seqs))
        if key not in cache:
            path = model_path(legacy_results, legacy_prefix, name)
            if not os.path.exists(path):
                pytest.skip("model missing: " + path)
            model = model_io.load_model(path)
            # the goldens are CPU float32; 'highest' keeps GPUs from using TF32 (traced inside the context)
            with jax.default_matmul_precision("highest"):
                logits, logcounts = model.predict(seqs.astype(np.float32), batch_size=32, verbose=0)
            cache[key] = (logits, logcounts)
        return cache[key]

    return run


def load_golden(goldens, name):
    path = os.path.join(goldens, "pred_{}.npz".format(name))
    if not os.path.exists(path):
        pytest.skip("golden missing: " + path)
    golden = np.load(path)
    seqs = np.load(os.path.join(goldens, "pred_inputs.npz"))["seqs"][:int(golden["n"])]
    return golden, seqs


@pytest.mark.parametrize("name", golden_names())
def test_predictions_match_legacy(name, goldens, predict, record_property):
    golden, seqs = load_golden(goldens, name)
    logits, logcounts = predict(name, seqs)
    assert logits.shape == golden["logits"].shape and logcounts.shape == golden["logcounts"].shape
    d_logits = float(np.abs(logits - golden["logits"]).max())
    d_logcounts = float(np.abs(logcounts - golden["logcounts"]).max())
    record_property("max_abs_diff_logits", d_logits)
    record_property("max_abs_diff_logcounts", d_logcounts)
    print("{} ({} seqs, {}): max|dlogits|={:.3e} max|dlogcounts|={:.3e}".format(
        name, len(seqs), jax.default_backend(), d_logits, d_logcounts))
    assert d_logits <= TOL and d_logcounts <= TOL, (name, d_logits, d_logcounts)


def test_full_model_is_nobias_plus_bias(goldens, predict, record_property):
    goldens_by_name = {name: load_golden(goldens, name) for name in FULL_MODELS}
    seqs = goldens_by_name["chrombpnet"][1]
    full_logits, full_logcounts = predict("chrombpnet", seqs)
    nb_logits, nb_logcounts = predict("chrombpnet_nobias", seqs)
    bias_logits, bias_logcounts = predict("bias_model_scaled", seqs)
    d_logits = float(np.abs(full_logits - (nb_logits + bias_logits)).max())
    combined = logsumexp(np.concatenate([nb_logcounts, bias_logcounts], axis=1), axis=1, keepdims=True)
    d_logcounts = float(np.abs(full_logcounts - combined).max())
    record_property("max_abs_diff_logits", d_logits)
    record_property("max_abs_diff_logcounts", d_logcounts)
    print("chrombpnet vs nobias+bias ({} seqs): max|dlogits|={:.3e} max|dlogcounts|={:.3e}".format(
        len(seqs), d_logits, d_logcounts))
    assert d_logits <= TOL and d_logcounts <= TOL, (d_logits, d_logcounts)
