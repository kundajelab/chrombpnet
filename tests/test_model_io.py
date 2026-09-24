"""model_io.load_model on hand-built TF-Keras 2.x (.h5, keras_version 2.x) files: nested-model node indices and the
1.x logsumexp Lambda count head, and re-saving such a model with Keras 3."""
import base64
import copy
import json
import marshal

import numpy as np
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)
import h5py
import keras
from keras.src.legacy.saving import legacy_h5_format

from chrombpnet.training.utils import model_io
from chrombpnet.training.utils.layers import LogSumExp, LogSumExpCompat

TF_KERAS_VERSION = "2.12.0"
N = 3


# ---- TF-Keras 2.x config builders (the JSON a 1.x .h5 file holds in its model_config attribute) ----

def input_cfg(name, n=N):
    return {"class_name": "InputLayer", "name": name, "inbound_nodes": [],
            "config": {"batch_input_shape": [None, n], "dtype": "float32", "sparse": False, "ragged": False,
                       "name": name}}


def node(*refs):
    """One legacy inbound node: the call's inputs as [layer, node_index, tensor_index, kwargs]."""
    return [[name, node_index, tensor_index, {}] for name, node_index, tensor_index in refs]


def dense_cfg(name, units, *nodes):
    return {"class_name": "Dense", "name": name, "inbound_nodes": list(nodes),
            "config": {"name": name, "trainable": True, "dtype": "float32", "units": units, "activation": "linear",
                       "use_bias": True, "kernel_initializer": {"class_name": "GlorotUniform", "config": {"seed": None}},
                       "bias_initializer": {"class_name": "Zeros", "config": {}}, "kernel_regularizer": None,
                       "bias_regularizer": None, "activity_regularizer": None, "kernel_constraint": None,
                       "bias_constraint": None}}


def merge_cfg(class_name, name, *nodes, **config):
    return {"class_name": class_name, "name": name, "inbound_nodes": list(nodes),
            "config": dict({"name": name, "trainable": True, "dtype": "float32"}, **config)}


def functional_cfg(name, layers, input_layers, output_layers, *nodes):
    return {"class_name": "Functional", "name": name, "inbound_nodes": list(nodes),
            "config": {"name": name, "trainable": True, "layers": layers, "input_layers": input_layers,
                       "output_layers": output_layers}}


def sequential_cfg(name, layers, *nodes):
    return {"class_name": "Sequential", "name": name, "inbound_nodes": list(nodes),
            "config": {"name": name, "layers": layers}}


def logsumexp_lambda_cfg(name, *nodes):
    """The 1.x `Lambda(lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True))` count head."""
    fn = lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)  # noqa: E731, F821 (never called)
    code = base64.encodebytes(marshal.dumps(fn.__code__)).decode("ascii")
    return {"class_name": "Lambda", "name": name, "inbound_nodes": list(nodes),
            "config": {"name": name, "trainable": True, "dtype": "float32", "function": [code, None, None],
                       "function_type": "lambda", "module": "chrombpnet.training.models.chrombpnet_with_bias_model",
                       "output_shape": None, "output_shape_type": "raw", "output_shape_module": None,
                       "arguments": {}}}


def write_legacy_h5(path, model_config, weights_from):
    """A TF-Keras 2.x style .h5: model_config JSON + keras_version 2.x, weights (by layer) from a Keras 3 model
    with the same layers."""
    with h5py.File(path, "w") as f:
        f.attrs["model_config"] = json.dumps(dict(model_config, keras_version=TF_KERAS_VERSION,
                                                  backend="tensorflow"))
        f.attrs["keras_version"] = TF_KERAS_VERSION
        f.attrs["backend"] = "tensorflow"
        legacy_h5_format.save_weights_to_hdf5_group(f.create_group("model_weights"), weights_from)
    return str(path)


def predict(model, x):
    out = model.predict(x, verbose=0)
    return [np.asarray(o) for o in (out if isinstance(out, list) else [out])]


def assert_same_outputs(got, want):
    assert len(got) == len(want)
    for g, w in zip(got, want):
        np.testing.assert_allclose(g, w, rtol=1e-6, atol=1e-6)


# ---- chrombpnet.h5-like: sequence -> model_wo_bias, model (each called once) -> Add / Concatenate + Lambda ----

def bpnet_like(name, input_name):
    inp = keras.Input((N,), name=input_name)
    return keras.Model(inp, [keras.layers.Dense(4, name=name + "_profile")(inp),
                             keras.layers.Dense(1, name=name + "_counts")(inp)], name=name)


def bpnet_like_cfg(name, input_name, *nodes):
    return functional_cfg(name, [input_cfg(input_name), dense_cfg(name + "_profile", 4, node((input_name, 0, 0))),
                                 dense_cfg(name + "_counts", 1, node((input_name, 0, 0)))],
                          [[input_name, 0, 0]], [[name + "_profile", 0, 0], [name + "_counts", 0, 0]], *nodes)


def chrombpnet_like_config():
    # TF-Keras numbers each nested model's first outer call node 1 (node 0 is its construction)
    return {"class_name": "Functional", "config": {"name": "chrombpnet", "trainable": True, "layers": [
        input_cfg("sequence"),
        bpnet_like_cfg("model_wo_bias", "sequence", node(("sequence", 0, 0))),
        bpnet_like_cfg("model", "bias_sequence", node(("sequence", 0, 0))),
        merge_cfg("Add", "logits_profile_predictions", node(("model_wo_bias", 1, 0), ("model", 1, 0))),
        merge_cfg("Concatenate", "concatenate", node(("model_wo_bias", 1, 1), ("model", 1, 1)), axis=-1),
        logsumexp_lambda_cfg("logcount_predictions", node(("concatenate", 0, 0)))],
        "input_layers": [["sequence", 0, 0]],
        "output_layers": [["logits_profile_predictions", 0, 0], ["logcount_predictions", 0, 0]]}}


def chrombpnet_like_reference():
    seq = keras.Input((N,), name="sequence")
    nobias, bias = bpnet_like("model_wo_bias", "sequence"), bpnet_like("model", "bias_sequence")
    (nb_profile, nb_counts), (b_profile, b_counts) = nobias(seq), bias(seq)
    profile = keras.layers.Add(name="logits_profile_predictions")([nb_profile, b_profile])
    counts = LogSumExp(name="logcount_predictions")(keras.layers.Concatenate(name="concatenate")([nb_counts,
                                                                                                   b_counts]))
    return keras.Model(seq, [profile, counts], name="chrombpnet")


def randomize(model, seed=0):
    rng = np.random.RandomState(seed)
    model.set_weights([rng.normal(size=w.shape).astype("float32") for w in model.get_weights()])
    return model


@pytest.fixture(scope="module")
def chrombpnet_like_h5(tmp_path_factory):
    ref = randomize(chrombpnet_like_reference())
    path = write_legacy_h5(tmp_path_factory.mktemp("legacy") / "chrombpnet.h5", chrombpnet_like_config(), ref)
    return path, ref


def test_legacy_chrombpnet_like_nested_single_call(chrombpnet_like_h5):
    path, ref = chrombpnet_like_h5
    model = model_io.load_model(path)
    assert model.output_names == ["logits_profile_predictions", "logcount_predictions"]
    assert type(model.get_layer("logcount_predictions")) is LogSumExp
    x = np.random.RandomState(1).normal(size=(5, N)).astype("float32")
    got = predict(model, x)
    assert_same_outputs(got, predict(ref, x))
    # the heads combine both nested models: profile = nobias + bias, counts = logsumexp(nobias, bias)
    nb, b = predict(model.get_layer("model_wo_bias"), x), predict(model.get_layer("model"), x)
    np.testing.assert_allclose(got[0], nb[0] + b[0], rtol=1e-6, atol=1e-6)
    np.testing.assert_allclose(got[1], np.logaddexp(nb[1], b[1]), rtol=1e-6, atol=1e-6)


def test_legacy_chrombpnet_like_resaves(chrombpnet_like_h5, tmp_path):
    path, ref = chrombpnet_like_h5
    x = np.random.RandomState(2).normal(size=(3, N)).astype("float32")
    model = model_io.load_model(path)
    for name in ("migrated.h5", "migrated.keras"):
        model.save(str(tmp_path / name))
        again = model_io.load_model(str(tmp_path / name))
        assert type(again.get_layer("logcount_predictions")) is LogSumExp, name
        assert_same_outputs(predict(again, x), predict(ref, x))


# ---- a nested model called twice, both calls referenced (outputs must not swap) ----

def multicall_config(outputs):
    return {"class_name": "Functional", "config": {"name": "outer", "trainable": True, "layers": [
        input_cfg("a"), input_cfg("b"),
        functional_cfg("sub", [input_cfg("sub_in"), dense_cfg("d", 2, node(("sub_in", 0, 0)))],
                       [["sub_in", 0, 0]], [["d", 0, 0]], node(("a", 0, 0)), node(("b", 0, 0))),
        # consumes the two calls in reverse order, through legacy inbound_nodes refs
        merge_cfg("Subtract", "diff", node(("sub", 2, 0), ("sub", 1, 0)))],
        "input_layers": [["a", 0, 0], ["b", 0, 0]], "output_layers": outputs}}


def multicall_reference():
    a, b = keras.Input((N,), name="a"), keras.Input((N,), name="b")
    sub_in = keras.Input((N,), name="sub_in")
    sub = keras.Model(sub_in, keras.layers.Dense(2, name="d")(sub_in), name="sub")
    sa, sb = sub(a), sub(b)
    return keras.Model([a, b], [sa, sb, keras.layers.Subtract(name="diff")([sb, sa])], name="outer")


def test_legacy_nested_model_called_twice(tmp_path):
    ref = randomize(multicall_reference())
    # TF-Keras numbers the two calls 1 and 2
    path = write_legacy_h5(tmp_path / "multicall.h5", multicall_config([["sub", 1, 0], ["sub", 2, 0], ["diff", 0, 0]]),
                           ref)
    model = model_io.load_model(path)
    xa, xb = np.ones((2, N), "float32"), -2 * np.ones((2, N), "float32")
    got = predict(model, [xa, xb])
    want = predict(ref, [xa, xb])
    assert not np.allclose(want[0], want[1])
    assert_same_outputs(got, want)


def test_shift_leaves_functional_output_refs_to_keras():
    config = multicall_config([["sub", 1, 0], ["sub", 2, 0]])
    assert model_io._shift_nested_node_indices(config)
    layers = {layer["name"]: layer for layer in config["config"]["layers"]}
    assert layers["diff"]["inbound_nodes"] == [node(("sub", 1, 0), ("sub", 0, 0))]
    # Keras 3 get_tensor subtracts 1 itself for Functional refs in input_layers / output_layers
    assert config["config"]["output_layers"] == [["sub", 1, 0], ["sub", 2, 0]]


# ---- nested Sequential models (Keras 3 get_tensor does not shift those) ----

def sequential_reference(with_input):
    a, b = keras.Input((N,), name="a"), keras.Input((N,), name="b")
    layers = [keras.layers.Dense(2, name="sd")]
    seq = keras.Sequential(([keras.Input((N,))] if with_input else []) + layers, name="seq")
    sa, sb = seq(a), seq(b)
    return keras.Model([a, b], [sa, sb, keras.layers.Dense(1, name="head")(sb)], name="outer")


def sequential_config(with_input):
    # built with an input shape, TF-Keras gave the Sequential a construction node 0 (calls are 1 and 2);
    # without one it did not (calls are 0 and 1)
    first = 1 if with_input else 0
    layers = ([input_cfg("seq_in")] if with_input else []) + [dense_cfg("sd", 2)]
    return {"class_name": "Functional", "config": {"name": "outer", "trainable": True, "layers": [
        input_cfg("a"), input_cfg("b"),
        sequential_cfg("seq", layers, node(("a", 0, 0)), node(("b", 0, 0))),
        dense_cfg("head", 1, node(("seq", first + 1, 0)))],
        "input_layers": [["a", 0, 0], ["b", 0, 0]],
        "output_layers": [["seq", first, 0], ["seq", first + 1, 0], ["head", 0, 0]]}}


@pytest.mark.parametrize("with_input", [True, False], ids=["with_input_shape", "without_input_shape"])
def test_legacy_nested_sequential(tmp_path, with_input):
    ref = randomize(sequential_reference(with_input))
    path = write_legacy_h5(tmp_path / "seq.h5", sequential_config(with_input), ref)
    model = model_io.load_model(path)
    xa, xb = np.ones((2, N), "float32"), -2 * np.ones((2, N), "float32")
    want = predict(ref, [xa, xb])
    assert not np.allclose(want[0], want[1])
    assert_same_outputs(predict(model, [xa, xb]), want)


def test_shift_sequential_refs():
    with_input = sequential_config(True)
    assert model_io._shift_nested_node_indices(with_input)
    assert with_input["config"]["output_layers"] == [["seq", 0, 0], ["seq", 1, 0], ["head", 0, 0]]
    assert with_input["config"]["layers"][3]["inbound_nodes"] == [node(("seq", 1, 0))]
    without_input = sequential_config(False)
    unchanged = copy.deepcopy(without_input)
    assert not model_io._shift_nested_node_indices(without_input)
    assert without_input == unchanged


# ---- no nested model: nothing to renumber ----

def test_legacy_without_nested_models(tmp_path):
    config = {"class_name": "Functional", "config": {"name": "plain", "trainable": True, "layers": [
        input_cfg("a"), dense_cfg("d1", 4, node(("a", 0, 0))), dense_cfg("d2", 2, node(("d1", 0, 0))),
        dense_cfg("d3", 1, node(("d1", 0, 0)))],
        "input_layers": [["a", 0, 0]], "output_layers": [["d2", 0, 0], ["d3", 0, 0]]}}
    unchanged = copy.deepcopy(config)
    assert not model_io._shift_nested_node_indices(config)
    assert config == unchanged

    a = keras.Input((N,), name="a")
    h = keras.layers.Dense(4, name="d1")(a)
    ref = randomize(keras.Model(a, [keras.layers.Dense(2, name="d2")(h), keras.layers.Dense(1, name="d3")(h)]))
    path = write_legacy_h5(tmp_path / "plain.h5", config, ref)
    assert not model_io._is_legacy_h5_with_nested_models(path)
    x = np.random.RandomState(3).normal(size=(4, N)).astype("float32")
    assert_same_outputs(predict(model_io.load_model(path), x), predict(ref, x))


# ---- the 1.x logsumexp Lambda on its own: load -> save -> load ----

def lambda_only_config():
    return {"class_name": "Functional", "config": {"name": "model", "trainable": True, "layers": [
        input_cfg("concatenate", 2), logsumexp_lambda_cfg("logcount_predictions", node(("concatenate", 0, 0)))],
        "input_layers": [["concatenate", 0, 0]], "output_layers": [["logcount_predictions", 0, 0]]}}


@pytest.mark.parametrize("resave_as", ["h5", "keras"])
def test_legacy_lambda_round_trip(tmp_path, resave_as):
    inp = keras.Input((2,), name="concatenate")
    ref = keras.Model(inp, LogSumExp(name="logcount_predictions")(inp), name="model")
    path = write_legacy_h5(tmp_path / "legacy_lambda.h5", lambda_only_config(), ref)
    x = np.array([[0.0, 0.0], [1.0, -3.0]], "float32")
    want = np.logaddexp(x[:, :1], x[:, 1:])

    loaded = model_io.load_model(path)
    assert type(loaded.get_layer("logcount_predictions")) is LogSumExp
    np.testing.assert_allclose(predict(loaded, x)[0], want, rtol=1e-6)
    migrated = str(tmp_path / ("migrated." + resave_as))
    loaded.save(migrated)
    again = model_io.load_model(migrated)
    assert type(again.get_layer("logcount_predictions")) is LogSumExp
    np.testing.assert_allclose(predict(again, x)[0], want, rtol=1e-6)
    again.save(str(tmp_path / ("twice." + resave_as)))
    assert type(model_io.load_model(str(tmp_path / ("twice." + resave_as))).get_layer(
        "logcount_predictions")) is LogSumExp


@pytest.mark.parametrize("suffix", ["h5", "keras"])
def test_files_resaved_with_logsumexpcompat_still_load(tmp_path, suffix):
    # what re-saving a loaded 1.x model wrote before LogSumExpCompat.from_config returned a LogSumExp
    inp = keras.Input((2,), name="concatenate")
    model = keras.Model(inp, LogSumExpCompat(name="logcount_predictions")(inp))
    path = str(tmp_path / ("old_migrated." + suffix))
    model.save(path)
    loaded = model_io.load_model(path)
    assert isinstance(loaded.get_layer("logcount_predictions"), LogSumExp)
    x = np.array([[0.5, 2.0]], "float32")
    np.testing.assert_allclose(predict(loaded, x)[0], np.logaddexp(x[:, :1], x[:, 1:]), rtol=1e-6)


def test_other_lambdas_are_refused():
    config = logsumexp_lambda_cfg("my_lambda")["config"]
    with pytest.raises(ValueError, match="my_lambda"):
        LogSumExpCompat.from_config(config)
