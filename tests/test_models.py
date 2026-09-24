"""Bias (bpnet_model) and ChromBPNet (chrombpnet_with_bias_model) architecture files on Keras 3 / JAX."""
import types

import numpy as np
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)
import keras

import chrombpnet.training.models.bpnet_model as bpnet_model
import chrombpnet.training.models.chrombpnet_with_bias_model as chrombpnet_with_bias_model
from chrombpnet.training import runtime
from chrombpnet.training.utils import model_io
from chrombpnet.training.utils.layers import LogSumExp

INPUTLEN, OUTPUTLEN, FILTERS, N_DIL = 514, 200, 8, 2
COUNTS_LOSS_WEIGHT = 10.0
LOG_KEYS = {"loss", "logits_profile_predictions_loss", "logcount_predictions_loss"}


def legacy_args(seed=1234, **kw):
    # the attributes chrombpnet 1.x passed; the new optimizer / precision options are optional
    return types.SimpleNamespace(seed=seed, learning_rate=1e-3, **kw)


def model_params(**kw):
    params = {"filters": str(FILTERS), "n_dil_layers": str(N_DIL), "counts_loss_weight": str(COUNTS_LOSS_WEIGHT),
              "inputlen": str(INPUTLEN), "outputlen": str(OUTPUTLEN)}
    params.update(kw)
    return params


def random_batch(n=8, seed=0):
    rng = np.random.RandomState(seed)
    x = np.eye(4, dtype=np.int8)[rng.randint(0, 4, (n, INPUTLEN))]
    y = rng.poisson(0.5, (n, OUTPUTLEN)).astype(np.float64)
    return x, (y, np.log(1 + y.sum(-1, keepdims=True)))


@pytest.fixture(scope="module")
def bias_h5(tmp_path_factory):
    path = tmp_path_factory.mktemp("bias") / "bias.h5"
    bpnet_model.getModelGivenModelOptionsAndWeightInits(legacy_args(seed=7), model_params()).save(str(path))
    return str(path)


@pytest.fixture(scope="module")
def chrombpnet_model(bias_h5):
    return chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
        legacy_args(), model_params(bias_model_path=bias_h5))


def bias_submodel(model):
    return [layer for layer in model.layers if isinstance(layer, keras.Model) and layer.name != "model_wo_bias"][0]


def test_bias_model_contract():
    model = bpnet_model.getModelGivenModelOptionsAndWeightInits(legacy_args(), model_params())
    assert model.output_names == ["logits_profile_predictions", "logcount_predictions"]
    assert [tuple(o.shape) for o in model.outputs] == [(None, OUTPUTLEN), (None, 1)]
    # find_chrombpnet_hyperparams shifts the bias of the last layer, which must be the count Dense
    assert model.layers[-1].name == "logcount_predictions" and isinstance(model.layers[-1], keras.layers.Dense)
    names = {layer.name for layer in model.layers}
    assert {"sequence", "bpnet_1st_conv", "bpnet_1conv", "bpnet_2conv", "bpnet_1crop", "bpnet_2crop",
            "prof_out_precrop", "logits_profile_predictions_preflatten", "gap"} <= names
    assert isinstance(model.optimizer, keras.optimizers.Adam)
    assert float(model.optimizer.learning_rate) == pytest.approx(1e-3)
    assert not model.optimizer.use_ema
    bpnet_model.save_model_without_bias(model, "unused")  # a no-op for the bias architecture


def test_seed_sets_initial_weights():
    w1 = bpnet_model.getModelGivenModelOptionsAndWeightInits(legacy_args(seed=3), model_params()).get_weights()
    w2 = bpnet_model.getModelGivenModelOptionsAndWeightInits(legacy_args(seed=3), model_params()).get_weights()
    w3 = bpnet_model.getModelGivenModelOptionsAndWeightInits(legacy_args(seed=4), model_params()).get_weights()
    assert all(np.array_equal(a, b) for a, b in zip(w1, w2))
    assert not all(np.array_equal(a, b) for a, b in zip(w1, w3))


def test_bias_fit_log_keys_and_roundtrip(tmp_path):
    model = bpnet_model.getModelGivenModelOptionsAndWeightInits(legacy_args(), model_params())
    x, y = random_batch()
    history = model.fit(x, y, batch_size=4, epochs=1, verbose=0, shuffle=False)
    logs = {k: v[-1] for k, v in history.history.items()}
    assert set(logs) == LOG_KEYS
    assert np.isfinite(list(logs.values())).all()
    # per-output losses are unweighted; the total is the weighted sum
    assert logs["loss"] == pytest.approx(
        logs["logits_profile_predictions_loss"] + COUNTS_LOSS_WEIGHT * logs["logcount_predictions_loss"], rel=1e-5)

    path = str(tmp_path / "bias.h5")
    model.save(path)
    loaded = model_io.load_model(path)
    assert loaded.output_names == model.output_names
    for a, b in zip(model.predict(x, verbose=0), loaded.predict(x, verbose=0)):
        np.testing.assert_array_equal(a, b)


def test_chrombpnet_contract(chrombpnet_model):
    model = chrombpnet_model
    assert model.output_names == ["logits_profile_predictions", "logcount_predictions"]
    assert isinstance(model.get_layer("logcount_predictions"), LogSumExp)
    assert not any(isinstance(layer, keras.layers.Lambda) for layer in model.layers)
    nobias = model.get_layer("model_wo_bias")
    assert nobias.output_names == ["wo_bias_bpnet_logits_profile_predictions", "wo_bias_bpnet_logcount_predictions"]
    # only the model without bias trains; the pretrained bias model is frozen
    bias = bias_submodel(model)
    assert len(model.trainable_variables) == len(nobias.trainable_variables) > 0
    assert all(v.path.startswith("wo_bias_") for v in model.trainable_variables)
    assert bias.trainable_variables == [] and len(bias.non_trainable_variables) == len(bias.weights) > 0


def test_chrombpnet_heads_combine_bias_and_nobias(chrombpnet_model):
    model = chrombpnet_model
    x, _ = random_batch(4, seed=1)
    logits, logcounts = model.predict(x, verbose=0)
    nb_logits, nb_logcounts = model.get_layer("model_wo_bias").predict(x, verbose=0)
    b_logits, b_logcounts = bias_submodel(model).predict(x, verbose=0)
    np.testing.assert_allclose(logits, nb_logits + b_logits, rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(logcounts, np.logaddexp(nb_logcounts, b_logcounts), rtol=1e-6, atol=1e-6)


def test_chrombpnet_seeding_before_nobias(bias_h5):
    # the no-bias initial weights follow --seed
    build = lambda seed: chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
        legacy_args(seed=seed), model_params(bias_model_path=bias_h5)).get_layer("model_wo_bias").get_weights()
    w1, w2, w3 = build(5), build(5), build(6)
    assert all(np.array_equal(a, b) for a, b in zip(w1, w2))
    assert not all(np.array_equal(a, b) for a, b in zip(w1, w3))


def test_chrombpnet_fit_save_and_nobias(bias_h5, tmp_path):
    model = chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
        legacy_args(), model_params(bias_model_path=bias_h5))
    bias_before = bias_submodel(model).get_weights()
    nobias_before = model.get_layer("model_wo_bias").get_weights()
    x, y = random_batch()
    history = model.fit(x, y, batch_size=4, epochs=1, verbose=0)
    assert set(history.history) == LOG_KEYS
    assert all(np.array_equal(a, b) for a, b in zip(bias_before, bias_submodel(model).get_weights()))
    assert not all(np.array_equal(a, b) for a, b in zip(nobias_before, model.get_layer("model_wo_bias").get_weights()))

    prefix = str(tmp_path / "chrombpnet")
    model.save(prefix + ".h5")
    chrombpnet_with_bias_model.save_model_without_bias(model, prefix)

    full = model_io.load_model(prefix + ".h5")
    assert full.output_names == model.output_names
    assert isinstance(full.get_layer("logcount_predictions"), LogSumExp)
    assert len(full.trainable_variables) == len(model.trainable_variables)
    for a, b in zip(model.predict(x, verbose=0), full.predict(x, verbose=0)):
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-6)

    nobias = model_io.load_model(prefix + "_nobias.h5")
    assert nobias.name == "model_wo_bias"
    assert nobias.output_names == ["wo_bias_bpnet_logits_profile_predictions", "wo_bias_bpnet_logcount_predictions"]
    for a, b in zip(model.get_layer("model_wo_bias").predict(x, verbose=0), nobias.predict(x, verbose=0)):
        np.testing.assert_allclose(a, b, rtol=1e-6, atol=1e-6)


def test_bf16_keeps_heads_and_bias_float32(bias_h5, tmp_path):
    previous = keras.config.dtype_policy()
    try:
        runtime.configure_precision("bf16")
        bias = bpnet_model.getModelGivenModelOptionsAndWeightInits(legacy_args(), model_params())
        assert [o.dtype for o in bias.outputs] == ["float32", "float32"]
        assert bias.get_layer("bpnet_1conv").compute_dtype == "bfloat16"
        for name in ("prof_out_precrop", "logits_profile_predictions_preflatten", "logits_profile_predictions", "gap",
                     "logcount_predictions"):
            assert bias.get_layer(name).compute_dtype == "float32", name

        model = chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
            legacy_args(), model_params(bias_model_path=bias_h5))
        assert [o.dtype for o in model.outputs] == ["float32", "float32"]
        nobias = model.get_layer("model_wo_bias")
        assert [o.dtype for o in nobias.outputs] == ["float32", "float32"]
        assert nobias.get_layer("wo_bias_bpnet_2conv").compute_dtype == "bfloat16"
        assert {layer.dtype_policy.name for layer in bias_submodel(model).layers} == {"float32"}
        assert all(v.dtype == "float32" for v in model.weights)
        x, y = random_batch(4)
        history = model.fit(x, y, batch_size=4, epochs=1, verbose=0)
        assert np.isfinite(history.history["loss"]).all()

        # what train.py does before saving a model trained in bf16
        runtime.set_float32_policy(model)
        chrombpnet_with_bias_model.save_model_without_bias(model, str(tmp_path / "m"))
    finally:
        keras.config.set_dtype_policy(previous)
    loaded = model_io.load_model(str(tmp_path / "m_nobias.h5"))
    assert {layer.dtype_policy.name for layer in loaded.layers} == {"float32"}


def test_float32_policy_context_restores_training_policies(bias_h5, tmp_path):
    # what the bf16 ModelCheckpoint does at each save: float32 file, then back to the per-layer bf16 policies
    previous = keras.config.dtype_policy()
    try:
        runtime.configure_precision("bf16")
        model = chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
            legacy_args(), model_params(bias_model_path=bias_h5))
    finally:
        keras.config.set_dtype_policy(previous)
    all_layers = lambda m: m._flatten_layers(include_self=True, recursive=True)
    before = [(layer.name, layer.dtype_policy.name) for layer in all_layers(model)]
    assert {"mixed_bfloat16", "float32"} <= {name for _, name in before}
    with runtime.float32_policy(model):
        assert {layer.dtype_policy.name for layer in all_layers(model)} == {"float32"}
        model.save(str(tmp_path / "checkpoint.h5"))
    assert [(layer.name, layer.dtype_policy.name) for layer in all_layers(model)] == before
    loaded = model_io.load_model(str(tmp_path / "checkpoint.h5"))
    assert {layer.dtype_policy.name for layer in all_layers(loaded)} == {"float32"}
