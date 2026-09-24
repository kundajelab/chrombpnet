"""Training parity with chrombpnet 1.x (TF-Keras), against goldens from tests/goldens/make_goldens.py.

CHROMBPNET_GOLDENS: the goldens directory (trace_<name>/{init.h5,final.h5,batches.npz,trace.json}).
CHROMBPNET_LEGACY_RESULTS / CHROMBPNET_LEGACY_PREFIX: a chrombpnet 1.x results tree and its run name, for
full_models/<prefix>/models/bias_model_scaled.h5 (see tests/conftest.py and tests/goldens/run_goldens.sh).
Each test skips when its variable is unset or the files are missing.
"""
import json
import os
import types

import numpy as np
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)
import jax
import keras

import chrombpnet.training.models.chrombpnet_with_bias_model as chrombpnet_with_bias_model
from chrombpnet.training.utils import model_io
from chrombpnet.training.utils.losses import multinomial_nll

LOG_KEYS = ("loss", "logits_profile_predictions_loss", "logcount_predictions_loss")


def env_path(var, *parts):
    root = os.environ.get(var)
    if not root:
        pytest.skip("{} is not set".format(var))
    path = os.path.join(root, *parts)
    if not os.path.exists(path):
        pytest.skip("missing {}".format(path))
    return path


@pytest.fixture
def highest_precision():
    # the goldens were made on CPU (no TF32); compare at full float32
    previous = jax.config.jax_default_matmul_precision
    jax.config.update("jax_default_matmul_precision", "highest")
    yield
    jax.config.update("jax_default_matmul_precision", previous)


def bias_submodel(model):
    return [layer for layer in model.layers if isinstance(layer, keras.Model) and layer.name != "model_wo_bias"][0]


def rel_l2(a, b):
    return np.linalg.norm(np.asarray(a, np.float64) - b) / max(np.linalg.norm(np.asarray(b, np.float64)), 1e-12)


@pytest.mark.parametrize("name", ["bias_128x4", "chrombpnet_32x4"])
def test_adam_trace_matches_legacy(name, highest_precision):
    tdir = env_path("CHROMBPNET_GOLDENS", "trace_" + name)
    meta = json.load(open(os.path.join(tdir, "trace.json")))
    params = meta["params"]
    batches = np.load(os.path.join(tdir, "batches.npz"))

    model = model_io.load_model(os.path.join(tdir, "init.h5"))
    # the output names make the log keys
    assert model.output_names == ["logits_profile_predictions", "logcount_predictions"]
    if params.get("bias_model_path"):
        # chrombpnet.h5 with the legacy logsumexp Lambda; the bias model must still be frozen after loading
        nobias = model.get_layer("model_wo_bias")
        assert len(model.trainable_variables) == len(nobias.trainable_variables), (
            "the legacy frozen bias model came back trainable from model_io.load_model")
    model.compile(optimizer=keras.optimizers.Adam(learning_rate=1e-3), loss=[multinomial_nll, "mse"],
                  loss_weights=[1, float(params["counts_loss_weight"])])

    bs = meta["batch_size"]
    x, y, logy = batches["x"].astype(np.float32), batches["y"], batches["logy"]
    steps = []
    for i in range(len(meta["steps"])):
        sl = slice(i * bs, (i + 1) * bs)
        # TF-Keras train_on_batch reset the metrics every call (per-batch values); Keras 3 accumulates
        model.reset_metrics()
        steps.append({k: float(v) for k, v in model.train_on_batch(x[sl], (y[sl], logy[sl]), return_dict=True).items()})

    for key in LOG_KEYS:
        assert key in meta["steps"][0] and key in steps[0], key
        legacy = np.array([s[key] for s in meta["steps"]])
        ours = np.array([s[key] for s in steps])
        rel = np.abs(ours - legacy) / np.abs(legacy)
        assert rel[0] < 1e-5, (key, ours[0], legacy[0])
        assert rel[-1] < 1e-3, (key, ours[-1], legacy[-1])

    final = model_io.load_model(os.path.join(tdir, "final.h5"))
    ours, legacy = model.get_weights(), final.get_weights()
    paths = [v.path for v in model.weights]
    assert [w.shape for w in ours] == [w.shape for w in legacy]
    everything = rel_l2(np.concatenate([w.ravel() for w in ours]), np.concatenate([w.ravel() for w in legacy]))
    assert everything < 5e-3, everything
    # The profile-head bias shifts every logit equally, so softmax (and the loss) ignore it: its true gradient is 0
    # and Adam's normalised steps on it follow float noise in either implementation. Exclude it.
    per_array = {p: rel_l2(a, b) for p, a, b in zip(paths, ours, legacy) if not p.endswith("prof_out_precrop/bias")}
    worst = max(per_array, key=per_array.get)
    assert per_array[worst] < 5e-2, (worst, per_array[worst])


def test_legacy_bias_as_pretrained_bias(tmp_path, legacy_prefix):
    bias_h5 = env_path("CHROMBPNET_LEGACY_RESULTS", "full_models", legacy_prefix, "models", "bias_model_scaled.h5")
    bias = model_io.load_model(bias_h5)
    inputlen, outputlen = bias.inputs[0].shape[1], bias.outputs[0].shape[1]
    params = {"filters": "8", "n_dil_layers": "2", "counts_loss_weight": "10", "inputlen": str(inputlen),
              "outputlen": str(outputlen), "bias_model_path": bias_h5}
    model = chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
        types.SimpleNamespace(seed=1234, learning_rate=1e-3), params)
    assert len(model.trainable_variables) == len(model.get_layer("model_wo_bias").trainable_variables)

    rng = np.random.RandomState(0)
    x = np.eye(4, dtype=np.int8)[rng.randint(0, 4, (4, inputlen))]
    y = rng.poisson(0.2, (4, outputlen)).astype(np.float32)
    model.fit(x, (y, np.log1p(y.sum(-1, keepdims=True))), batch_size=4, epochs=1, verbose=0)
    # still the legacy bias model after a training step
    for a, b in zip(bias_submodel(model).predict(x, verbose=0), bias.predict(x, verbose=0)):
        np.testing.assert_allclose(a, b, rtol=1e-5, atol=1e-5)

    chrombpnet_with_bias_model.save_model_without_bias(model, str(tmp_path / "m"))
    assert model_io.load_model(str(tmp_path / "m_nobias.h5")).name == "model_wo_bias"
