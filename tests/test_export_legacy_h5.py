"""`chrombpnet export --legacy-h5`: TF-Keras 2.x full-model .h5 files from Keras 3 bias / chrombpnet / no-bias models.

Checks the file layout against what TF-Keras 2.12 wrote for chrombpnet 1.x (the trace_* reference files in
CHROMBPNET_GOLDENS, when set), the round trip through model_io.load_model, the dataset reads of bpnet-lite's
BPNet.from_chrombpnet, and re-exporting the 1.x reference files. CHROMBPNET_TF_PYTHON (a Python interpreter with
TensorFlow 2.x; run with TF 2.8 / Python 3.9 and TF 2.12 / Python 3.10) additionally loads the exports with
TF-Keras itself, the way chrombpnet 1.x and the variant-scorer load models.
"""
import base64
import json
import os
import re
import subprocess
import sys
import textwrap
import types

import numpy as np
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)
import h5py
import keras

import chrombpnet.parsers as parsers
import chrombpnet.training.models.bpnet_model as bpnet_model
import chrombpnet.training.models.chrombpnet_with_bias_model as chrombpnet_with_bias_model
from chrombpnet.helpers.postprocessing.export_legacy_h5 import (LEGACY_KERAS_VERSION, LEGACY_LOGSUMEXP_BYTECODE,
                                                                 LEGACY_LOGSUMEXP_MODULE, export_legacy_h5,
                                                                 logsumexp_lambda_fields)
from chrombpnet.training.utils import model_io
from chrombpnet.training.utils.layers import LOGSUMEXP_LAMBDA_FUNCTION, LogSumExp, LogSumExpCompat

# the geometry of the chrombpnet 1.x reference files (legacy defaults), with few filters. The bias and no-bias
# models differ in width, so the weights of one cannot pass for the other's (the full model's outputs are symmetric
# in them: profile logits add up, counts go through logsumexp)
INPUTLEN, OUTPUTLEN, N_DIL = 2114, 1000, 4
FILTERS = {"bias": 16, "nobias": 8}
# filters of trace_bias_128x4 (and of the bias model of trace_chrombpnet_32x4) / of its no-bias model -> ours
GOLDEN_FILTERS = {128: FILTERS["bias"], 32: FILTERS["nobias"]}
BIAS_LAYER_NAMES = ["sequence", "bpnet_1st_conv", "bpnet_1conv", "bpnet_1crop", "add", "bpnet_2conv", "bpnet_2crop",
                    "add_1", "bpnet_3conv", "bpnet_3crop", "add_2", "bpnet_4conv", "bpnet_4crop", "add_3",
                    "prof_out_precrop", "logits_profile_predictions_preflatten", "gap", "logits_profile_predictions",
                    "logcount_predictions"]
NOBIAS_LAYER_NAMES = ["sequence", "wo_bias_bpnet_1st_conv", "wo_bias_bpnet_1conv", "wo_bias_bpnet_1crop", "add",
                      "wo_bias_bpnet_2conv", "wo_bias_bpnet_2crop", "add_1", "wo_bias_bpnet_3conv",
                      "wo_bias_bpnet_3crop", "add_2", "wo_bias_bpnet_4conv", "wo_bias_bpnet_4crop", "add_3",
                      "wo_bias_bpnet_prof_out_precrop", "wo_bias_bpnet_logitt_before_flatten", "gap",
                      "wo_bias_bpnet_logits_profile_predictions", "wo_bias_bpnet_logcount_predictions"]
FULL_LAYER_NAMES = ["sequence", "model_wo_bias", "model", "concatenate", "logits_profile_predictions",
                    "logcount_predictions"]
CONV_NAMES = ["bpnet_1st_conv"] + ["bpnet_{}conv".format(i) for i in range(1, N_DIL + 1)] + ["prof_out_precrop"]
BIAS_WEIGHT_LAYERS = CONV_NAMES + ["logcount_predictions"]
NOBIAS_WEIGHT_LAYERS = ["wo_bias_" + n for n in CONV_NAMES[:-1]] + ["wo_bias_bpnet_prof_out_precrop",
                                                                     "wo_bias_bpnet_logcount_predictions"]
TRACES = ["trace_bias_128x4", "trace_chrombpnet_32x4"]
# the count-head Lambda config of the 1.x chrombpnet.h5 files (the export's default count head)
LAMBDA_1X = {"name": "logcount_predictions", "trainable": True, "dtype": "float32",
             "function": [LEGACY_LOGSUMEXP_BYTECODE, None, None], "function_type": "lambda",
             "module": "chrombpnet.training.models.chrombpnet_with_bias_model", "output_shape": None,
             "output_shape_type": "raw", "output_shape_module": None, "arguments": {}}
LAMBDA_NAMED = dict(LAMBDA_1X, function=LOGSUMEXP_LAMBDA_FUNCTION, function_type="function", module=None)


def model_params(kind, **kw):
    params = {"filters": str(FILTERS[kind]), "n_dil_layers": str(N_DIL), "counts_loss_weight": "10.0",
              "inputlen": str(INPUTLEN), "outputlen": str(OUTPUTLEN)}
    params.update(kw)
    return params


def randomize(model, seed):
    rng = np.random.RandomState(seed)
    model.set_weights([(0.1 * rng.normal(size=w.shape)).astype("float32") for w in model.get_weights()])
    return model


def one_hot(n=3, seed=0):
    return np.eye(4, dtype="float32")[np.random.RandomState(seed).randint(0, 4, (n, INPUTLEN))]


def predict(model, x):
    return [np.asarray(o) for o in model.predict(x, verbose=0)]


# ---- models and exports ----

@pytest.fixture(scope="module")
def models(tmp_path_factory):
    """Keras 3 bias, chrombpnet and no-bias models built with the architecture files, and their Keras 3 files.
    The no-bias model is built after the bias model in the same process, so Keras 3 numbers its Add layers
    add_4 ... add_7 (TF-Keras 1.x files have add ... add_3). Bias: 16 filters, no-bias: 8 filters."""
    d = tmp_path_factory.mktemp("export_models")
    args = types.SimpleNamespace(seed=11, learning_rate=1e-3)
    bias = randomize(bpnet_model.getModelGivenModelOptionsAndWeightInits(args, model_params("bias")), 1)
    bias_h5 = str(d / "bias.h5")
    bias.save(bias_h5)
    full = chrombpnet_with_bias_model.getModelGivenModelOptionsAndWeightInits(
        args, model_params("nobias", bias_model_path=bias_h5))
    randomize(full.get_layer("model_wo_bias"), 2)
    full_h5 = str(d / "chrombpnet.h5")
    full.save(full_h5)
    chrombpnet_with_bias_model.save_model_without_bias(full, str(d / "chrombpnet"))
    return {"bias": bias, "full": full, "nobias": full.get_layer("model_wo_bias"),
            "paths": {"bias": bias_h5, "full": full_h5, "nobias": str(d / "chrombpnet_nobias.h5")}}


@pytest.fixture(scope="module")
def exported(models, tmp_path_factory):
    """Legacy exports: bias and no-bias from their Keras 3 .h5 files, the full model from the model object (with
    the default 1.x bytecode count head, and with the named-function one)."""
    d = tmp_path_factory.mktemp("exported")
    return {"bias": export_legacy_h5(models["paths"]["bias"], d / "bias.legacy.h5"),
            "nobias": export_legacy_h5(models["paths"]["nobias"], d / "nobias.legacy.h5"),
            "full": export_legacy_h5(models["full"], d / "chrombpnet.legacy.h5"),
            "full_named": export_legacy_h5(models["full"], d / "chrombpnet.named.legacy.h5", count_head="named")}


def trace_file(goldens_dir, trace, which="final.h5"):
    path = goldens_dir / trace / which
    if not path.exists():
        pytest.skip("{} not found".format(path))
    return str(path)


# ---- HDF5 helpers ----

def attr_kind(obj, key):
    """HDF5 type of an attribute: ('str', variable-length?, charset) / ('float', shape) / ..."""
    aid = h5py.h5a.open(obj.id, key.encode())
    t = aid.get_type()
    if isinstance(t, h5py.h5t.TypeStringID):
        return ("str", t.is_variable_str(), t.get_cset(), aid.get_space().shape)
    return (type(t).__name__, aid.get_space().shape)


def attr_values(obj):
    out = {}
    for key, value in obj.attrs.items():
        out[key] = [v.decode() if isinstance(v, bytes) else v for v in value.tolist()] \
            if isinstance(value, np.ndarray) else value
    return out


def weights_tree(path):
    """{path under model_weights: ('group', attrs, attr kinds) | ('dataset', shape, dtype)}."""
    tree = {}
    with h5py.File(path, "r") as f:
        group = f["model_weights"]
        tree[""] = ("group", attr_values(group), {k: attr_kind(group, k) for k in group.attrs})

        def visit(name, obj):
            if isinstance(obj, h5py.Dataset):
                tree[name] = ("dataset", obj.shape, obj.dtype.str)
            else:
                tree[name] = ("group", attr_values(obj), {k: attr_kind(obj, k) for k in obj.attrs})
        group.visititems(visit)
    return tree


def root_attrs(path):
    with h5py.File(path, "r") as f:
        return {k: attr_kind(f, k) for k in f.attrs}, list(f.keys())


def model_config(path):
    with h5py.File(path, "r") as f:
        return json.loads(f.attrs["model_config"])


def dataset_values(path):
    with h5py.File(path, "r") as f:
        out = {}
        f["model_weights"].visititems(
            lambda name, obj: out.__setitem__(name, obj[()]) if isinstance(obj, h5py.Dataset) else None)
        return out


def config_diff(a, b, path=""):
    """[(path, a, b)] where two JSON configs differ (dict key order included)."""
    if isinstance(a, dict) and isinstance(b, dict):
        out = [] if list(a) == list(b) else [(path + "{keys}", list(a), list(b))]
        return out + [d for k in a if k in b for d in config_diff(a[k], b[k], path + "/" + k)]
    if isinstance(a, list) and isinstance(b, list):
        out = [] if len(a) == len(b) else [(path + "{len}", len(a), len(b))]
        return out + [d for i, (x, y) in enumerate(zip(a, b)) for d in config_diff(x, y, "{}[{}]".format(path, i))]
    return [] if (a == b and type(a) is type(b)) else [(path, a, b)]


def walk_configs(config):
    """Every dict inside a JSON config."""
    if isinstance(config, dict):
        yield config
        for v in config.values():
            yield from walk_configs(v)
    elif isinstance(config, list):
        for v in config:
            yield from walk_configs(v)


def refs(config):
    """Every [layer, node_index, tensor_index, kwargs] inbound ref of a (possibly nested) model config."""
    for layer in config["config"]["layers"]:
        for node in layer["inbound_nodes"]:
            yield from node
        if layer["class_name"] == "Functional":
            yield from refs(layer)


def lambda_layer(config):
    return [layer for layer in config["config"]["layers"] if layer["class_name"] == "Lambda"]


def source_weights(kind, models):
    """{"<top layer>/<layer>/<weight>:0": value} of the Keras 3 source model, at the path TF-Keras gives each
    weight: in a full model the bias model's weights under `model`, the no-bias model's under `model_wo_bias`
    (independent of the exporter's naming code)."""
    parts = [("model", models["bias"]), ("model_wo_bias", models["nobias"])] if kind == "full" \
        else [(None, models[kind])]
    out = {}
    for top, model in parts:
        for layer in model.layers:
            for v in layer.weights:
                out["{}/{}/{}:0".format(top or layer.name, layer.name, v.name)] = keras.ops.convert_to_numpy(v)
    return out


def assert_variables_match(model, path, reference):
    """Every variable of a loaded model (of its top-level layers and of the layers of its nested models) equals
    the dataset model_weights/<top layer>/<layer>/<weight>:0 of `path` and the source weight of that path; every
    dataset belongs to a variable."""
    values = dataset_values(path)
    seen = []
    for top in model.layers:
        for sub in (top.layers if isinstance(top, keras.Model) else [top]):
            for v in sub.weights:
                key = "{}/{}/{}:0".format(top.name, sub.name, v.name)
                value = keras.ops.convert_to_numpy(v)
                np.testing.assert_array_equal(value, values[key], err_msg=key)
                np.testing.assert_array_equal(value, reference[key], err_msg=key)
                seen.append(key)
    assert sorted(seen) == sorted(values) == sorted(reference)


# ---- bpnet-lite ----

def bpnetlite_from_chrombpnet_reads(filename):
    """The h5py reads of bpnet-lite's BPNet.from_chrombpnet (jmschrei/bpnet-lite, bpnetlite/bpnet.py), copied
    without the torch conversion. Returns (n_layers, n_filters, {dataset path: array}).

    (from_chrombpnet_lite, with its `model_1` / conv1d_N naming, reads chrombpnet-lite files, a different format.
    ChromBPNet.from_chrombpnet calls BPNet.from_chrombpnet on the bias and the no-bias file.)"""
    reads = {}
    with h5py.File(filename, "r") as h5:
        w = h5['model_weights']

        if 'bpnet_1conv' in w.keys():
            prefix = ""
        else:
            prefix = "wo_bias_"

        def namer(prefix, suffix):
            return '{0}{1}/{0}{1}'.format(prefix, suffix)
        k, b = 'kernel:0', 'bias:0'

        def read(name):
            reads[name + "/" + k] = w[name][k][:]
            reads[name + "/" + b] = w[name][b][:]

        n_layers = 0
        for layer_name in w.keys():
            try:
                idx = int(layer_name.split("_")[-1].replace("conv", ""))
                n_layers = max(n_layers, idx)
            except ValueError:   # bpnet-lite: bare except around the same int() call
                pass

        name = namer(prefix, "bpnet_1conv")
        n_filters = w[name][k].shape[2]

        read(namer(prefix, 'bpnet_1st_conv'))
        for i in range(1, n_layers + 1):
            read(namer(prefix, 'bpnet_{}conv'.format(i)))

        prefix = prefix + "bpnet_" if prefix != "" else ""
        read(namer(prefix, 'prof_out_precrop'))
        read(namer(prefix, "logcount_predictions"))
    return n_layers, n_filters, reads


# ---- (1) structure ----

@pytest.mark.parametrize("kind", ["bias", "nobias", "full"])
def test_layout(kind, models, exported):
    path, model = exported[kind], models[kind]
    attrs, root_keys = root_attrs(path)
    # root: keras_version / backend as variable-length UTF-8 str, model_config as variable-length ASCII (what
    # TF-Keras 2.x save_model_to_hdf5 writes); no training_config / optimizer_weights (compile=False use)
    assert attrs == {"keras_version": ("str", True, h5py.h5t.CSET_UTF8, ()),
                     "backend": ("str", True, h5py.h5t.CSET_UTF8, ()),
                     "model_config": ("str", True, h5py.h5t.CSET_ASCII, ())}
    assert root_keys == ["model_weights"]
    with h5py.File(path, "r") as f:
        assert f.attrs["keras_version"] == LEGACY_KERAS_VERSION and f.attrs["backend"] == "tensorflow"

    tree = weights_tree(path)
    _, top_attrs, top_kinds = tree[""]
    layer_names = {"bias": BIAS_LAYER_NAMES, "nobias": NOBIAS_LAYER_NAMES, "full": FULL_LAYER_NAMES}[kind]
    assert top_attrs == {"layer_names": layer_names, "backend": "tensorflow", "keras_version": LEGACY_KERAS_VERSION}
    assert top_kinds["layer_names"] == ("str", True, h5py.h5t.CSET_ASCII, (len(layer_names),))
    groups = {n for n, v in tree.items() if v[0] == "group" and n and "/" not in n}
    assert groups == set(layer_names) | {"top_level_model_weights"}

    expected = {}   # group -> weight_names, in TF order
    if kind == "full":
        expected["model_wo_bias"] = [n + "/" + w for n in NOBIAS_WEIGHT_LAYERS for w in ("kernel:0", "bias:0")]
        expected["model"] = [n + "/" + w for n in BIAS_WEIGHT_LAYERS for w in ("kernel:0", "bias:0")]
    else:
        for n in (BIAS_WEIGHT_LAYERS if kind == "bias" else NOBIAS_WEIGHT_LAYERS):
            expected[n] = [n + "/kernel:0", n + "/bias:0"]
    for group in groups:
        weight_names = tree[group][1]["weight_names"]
        assert weight_names == expected.get(group, []), group
        if weight_names:
            assert tree[group][2]["weight_names"] == ("str", True, h5py.h5t.CSET_ASCII, (len(weight_names),))
        else:   # TF-Keras writes an empty float64 array
            assert tree[group][2]["weight_names"] == ("TypeFloatID", (0,))
    datasets = {n: v for n, v in tree.items() if v[0] == "dataset"}
    assert set(datasets) == {g + "/" + w for g, names in expected.items() for w in names}
    for name, (_, shape, dtype) in datasets.items():
        assert re.fullmatch(r"[^/]+/[^/]+/(kernel|bias):0", name), name
        assert dtype == "<f4"

    # every weight of the Keras model, by name
    values = dataset_values(path)
    leaves = [layer for layer in model._flatten_layers(include_self=False)
              if layer.weights and not isinstance(layer, keras.Model)]
    layers = {layer.name: layer for layer in leaves}
    assert len(layers) == len(leaves)
    for name, value in values.items():
        group, layer, weight = name.split("/")
        variable = {"kernel:0": 0, "bias:0": 1}[weight]
        np.testing.assert_array_equal(value, layers[layer].get_weights()[variable])
    assert len(values) == sum(len(layer.weights) for layer in layers.values())


@pytest.mark.parametrize("kind", ["bias", "nobias", "full"])
def test_model_config_is_keras2(kind, exported):
    config = model_config(exported[kind])
    assert set(config) == {"class_name", "config"} and config["class_name"] == "Functional"
    assert config["config"]["name"] == ("model_wo_bias" if kind == "nobias" else "model")
    for d in walk_configs(config):
        keras3_only = {"module", "registered_name", "build_config", "shared_object_id"} & set(d)
        if "function_type" in d:   # a Keras 2 Lambda config has a "module" key of its own
            keras3_only.discard("module")
        assert not keras3_only, d
        if "dtype" in d:
            assert d["dtype"] == "float32", d
        if d.get("class_name") == "InputLayer":
            assert list(d["config"]) == ["batch_input_shape", "dtype", "sparse", "ragged", "name"]
            assert d["config"]["batch_input_shape"] == [None, INPUTLEN, 4]
        if "inbound_nodes" in d:
            assert all(isinstance(node, list) for node in d["inbound_nodes"])   # list form, not Keras 3 dicts
    if kind == "full":
        # a nested model called once is node 1 (TF-Keras counts its construction as node 0)
        nested = {"model", "model_wo_bias"}
        top_refs = [ref for layer in config["config"]["layers"] for node in layer["inbound_nodes"] for ref in node]
        assert {ref[1] for ref in top_refs if ref[0] in nested} == {1}
        assert {ref[1] for ref in top_refs if ref[0] not in nested} == {0}
        bias = [layer for layer in config["config"]["layers"] if layer["name"] == "model"][0]
        # the frozen bias model: layers trainable=False inside a trainable nested model, as chrombpnet 1.x wrote
        assert bias["config"]["trainable"] is True
        assert {layer["config"].get("trainable") for layer in bias["config"]["layers"][1:]} == {False}
        # the count head: the Lambda of the 1.x files (Python 3.8 bytecode), keys in TF-Keras' order
        [lse] = lambda_layer(config)
        assert lse["name"] == "logcount_predictions"
        assert lse["config"] == LAMBDA_1X and list(lse["config"]) == list(LAMBDA_1X)
        assert lse["inbound_nodes"] == [[["concatenate", 0, 0, {}]]]
    else:
        assert not lambda_layer(config)
        assert config["config"]["layers"][-1]["class_name"] == "Dense"
        assert {ref[1] for ref in refs(config)} == {0}


def test_auto_names_become_tf_names(models, exported):
    # Keras 3 names: nested bias model 'functional', no-bias Add layers numbered after the bias model's
    source = models["full"]
    assert re.fullmatch(r"functional(_\d+)?", source.name)
    source_bias = [layer for layer in source.layers if isinstance(layer, keras.Model) and layer.name != "model_wo_bias"]
    assert re.fullmatch(r"functional(_\d+)?", source_bias[0].name)
    source_adds = [layer.name for layer in models["nobias"].layers if isinstance(layer, keras.layers.Add)]
    assert source_adds != ["add", "add_1", "add_2", "add_3"]
    config = model_config(exported["full"])
    assert [layer["name"] for layer in config["config"]["layers"]] == FULL_LAYER_NAMES
    for nested in config["config"]["layers"][1:3]:
        assert [layer["name"] for layer in nested["config"]["layers"] if layer["class_name"] == "Add"] == \
            ["add", "add_1", "add_2", "add_3"]


def test_named_count_head(models, exported):
    """count_head="named": the same file except for the count-head Lambda, which names its function."""
    named, default = model_config(exported["full_named"]), model_config(exported["full"])
    [lse] = lambda_layer(named)
    assert lse["config"] == LAMBDA_NAMED and list(lse["config"]) == list(LAMBDA_NAMED)
    assert config_diff(named, default) == [
        ("/config/layers[5]/config/" + k, LAMBDA_NAMED[k], LAMBDA_1X[k]) for k in ("function", "function_type", "module")]
    assert weights_tree(exported["full_named"]) == weights_tree(exported["full"])
    assert root_attrs(exported["full_named"]) == root_attrs(exported["full"])
    loaded = model_io.load_model(exported["full_named"])
    assert type(loaded.get_layer("logcount_predictions")) is LogSumExp
    assert_variables_match(loaded, exported["full_named"], source_weights("full", models))
    x = one_hot()
    for got, want in zip(predict(loaded, x), predict(models["full"], x)):
        np.testing.assert_allclose(got, want, rtol=0, atol=1e-6)


def test_bad_count_head_is_refused(models, tmp_path):
    with pytest.raises(ValueError, match="count_head must be one of bytecode, named"):
        export_legacy_h5(models["paths"]["full"], tmp_path / "bad.h5", count_head="source")
    assert not os.listdir(tmp_path)


def test_count_head_constant_is_the_1x_bytecode():
    # base64 of marshalled Python 3.8 code for `lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)`,
    # written from chrombpnet_with_bias_model.py (the reference files below check it byte for byte)
    code = base64.decodebytes(LEGACY_LOGSUMEXP_BYTECODE.encode("ascii"))
    assert base64.encodebytes(code).decode("ascii") == LEGACY_LOGSUMEXP_BYTECODE
    for name in (b"tf", b"math", b"reduce_logsumexp", b"axis", b"keepdims", b"chrombpnet_with_bias_model.py"):
        assert name in code
    assert LEGACY_LOGSUMEXP_MODULE == chrombpnet_with_bias_model.__name__
    assert logsumexp_lambda_fields() == {k: v for k, v in LAMBDA_1X.items() if k not in ("name", "trainable", "dtype")}


@pytest.mark.parametrize("which", ["init.h5", "final.h5"])
def test_count_head_is_the_1x_lambda(which, goldens_dir, exported):
    """The count-head Lambda of a Keras 3 model's export is the one of the 1.x chrombpnet.h5 files, read from
    them: config (keys in the same order) and inbound nodes."""
    golden = trace_file(goldens_dir, "trace_chrombpnet_32x4", which)
    [want] = lambda_layer(model_config(golden))
    [got] = lambda_layer(model_config(exported["full"]))
    assert json.dumps(got) == json.dumps(want)
    assert want["config"]["function"][0] == LEGACY_LOGSUMEXP_BYTECODE


def test_names_kept_without_normalization_or_on_clash(tmp_path):
    inp = keras.Input((3,), name="x")
    a, b = keras.layers.Dense(2, name="add")(inp), keras.layers.Dense(2, name="d2")(inp)
    add = keras.layers.Add(name="add_7")([a, b])
    model = keras.Model(inp, [add, keras.layers.Dense(1, name="out")(add)], name="functional_3")
    with pytest.warns(UserWarning, match="Keeping the Keras 3 layer names"):
        config = model_config(export_legacy_h5(model, tmp_path / "clash.h5"))
    # 'add_7' -> 'add' would clash with the Dense named 'add': names stay; the model itself is still renamed
    assert [layer["name"] for layer in config["config"]["layers"]] == ["x", "add", "d2", "add_7", "out"]
    assert config["config"]["name"] == "model"
    config = model_config(export_legacy_h5(model, tmp_path / "keep.h5", normalize_names=False))
    assert config["config"]["name"] == "functional_3"


def test_mixed_precision_layers_are_written_float32(tmp_path):
    inp = keras.Input((5,), name="x")
    out = keras.layers.Dense(2, name="d", dtype="mixed_bfloat16")(inp)
    out = keras.layers.Dense(2, name="b", dtype="bfloat16")(out)   # bfloat16 variables
    model = keras.Model(inp, keras.layers.Dense(1, name="out", dtype="float32")(out))
    path = export_legacy_h5(model, tmp_path / "bf16.h5")
    assert [layer["config"]["dtype"] for layer in model_config(path)["config"]["layers"]] == ["float32"] * 4
    assert {v.dtype for v in dataset_values(path).values()} == {np.dtype("float32")}


def test_unsupported_layers_are_refused(tmp_path):
    inp = keras.Input((5,), name="x")
    model = keras.Model(inp, keras.layers.Dense(1)(keras.layers.Dropout(0.1)(inp)))
    with pytest.raises(ValueError, match="not supported by the legacy .h5 export"):
        export_legacy_h5(model, tmp_path / "dropout.h5")
    assert not os.path.exists(tmp_path / "dropout.h5") and not os.path.exists(str(tmp_path / "dropout.h5") + ".tmp")


@pytest.mark.parametrize("trace,kind", [("trace_bias_128x4", "bias"), ("trace_chrombpnet_32x4", "full")])
def test_layout_matches_the_1x_reference_files(trace, kind, goldens_dir, exported):
    """A Keras 3 model's export has the layout of the TF-Keras 2.12 file of the same architecture: same
    attributes (and HDF5 types), layer_names, weight_names, groups, dataset paths and dtypes, and the same
    model_config (count-head Lambda included) except that shapes and `filters` differ by the number of filters."""
    golden, path = trace_file(goldens_dir, trace), exported[kind]
    golden_attrs, golden_keys = root_attrs(golden)
    attrs, keys = root_attrs(path)
    assert attrs == {k: v for k, v in golden_attrs.items() if k != "training_config"}
    assert keys == [k for k in golden_keys if k != "optimizer_weights"]

    want, got = weights_tree(golden), weights_tree(path)
    assert set(got) == set(want)
    for name in want:
        if want[name][0] == "group":
            assert got[name][1:] == want[name][1:], name
        else:
            shape = tuple(GOLDEN_FILTERS.get(s, s) for s in want[name][1])
            assert got[name] == ("dataset", shape, want[name][2]), name

    diffs = config_diff(model_config(path), model_config(golden))
    assert diffs and all(d[0].endswith("/filters") and GOLDEN_FILTERS[d[2]] == d[1] for d in diffs), diffs


# ---- (2) round trip ----

@pytest.mark.parametrize("kind", ["bias", "nobias", "full"])
def test_round_trip(kind, models, exported):
    x = one_hot()
    loaded = model_io.load_model(exported[kind])
    # each variable by path, against the file and the source model (predictions alone cannot tell a full model
    # from one with the bias and no-bias weights swapped)
    assert_variables_match(loaded, exported[kind], source_weights(kind, models))
    for got, want in zip(predict(loaded, x), predict(models[kind], x)):
        np.testing.assert_allclose(got, want, rtol=0, atol=1e-6)
    assert loaded.output_names == models[kind].output_names
    if kind == "full":
        assert type(loaded.get_layer("logcount_predictions")) is LogSumExp
        assert [layer.name for layer in loaded.layers] == FULL_LAYER_NAMES
        for nested, source in (("model_wo_bias", "nobias"), ("model", "bias")):
            for got, want in zip(predict(loaded.get_layer(nested), x), predict(models[source], x)):
                np.testing.assert_allclose(got, want, rtol=0, atol=1e-6, err_msg=nested)
        assert not loaded.get_layer("model").trainable_weights   # still frozen
        assert len(loaded.trainable_weights) == len(models["full"].trainable_weights)


def test_nested_model_object_exports_like_its_file(models, tmp_path):
    direct = export_legacy_h5(models["nobias"], tmp_path / "direct.h5")
    from_file = export_legacy_h5(models["paths"]["nobias"], tmp_path / "file.h5")
    assert weights_tree(direct) == weights_tree(from_file)
    assert model_config(direct) == model_config(from_file)


def test_keras_file_exports_like_the_model(models, exported, tmp_path):
    models["full"].save(str(tmp_path / "chrombpnet.keras"))
    path = export_legacy_h5(tmp_path / "chrombpnet.keras", tmp_path / "from_keras.h5")
    assert weights_tree(path) == weights_tree(exported["full"])
    assert model_config(path) == model_config(exported["full"])
    got, want = dataset_values(path), dataset_values(exported["full"])
    for name in want:
        np.testing.assert_array_equal(got[name], want[name])


def test_logsumexp_compat_accepts_both_count_heads():
    base = {"name": "logcount_predictions", "trainable": True, "dtype": "float32", "output_shape": None,
            "arguments": {}}
    # as Keras 3's legacy loader passes them (function_type / module dropped) and as written in the file
    for config in (dict(base, function=LOGSUMEXP_LAMBDA_FUNCTION), LAMBDA_NAMED,
                   dict(base, function=[LEGACY_LOGSUMEXP_BYTECODE, None, None]), LAMBDA_1X):
        assert type(LogSumExpCompat.from_config(config)) is LogSumExp
    for config in (dict(base, function="some_other_function", function_type="function"),
                   dict(base, function=LOGSUMEXP_LAMBDA_FUNCTION, name="other_head")):
        with pytest.raises(ValueError, match="Cannot load Lambda layer"):
            LogSumExpCompat.from_config(config)


# ---- (3) bpnet-lite ----

@pytest.mark.parametrize("kind", ["bias", "nobias"])
def test_bpnetlite_finds_every_weight(kind, models, exported):
    n_layers, n_filters, reads = bpnetlite_from_chrombpnet_reads(exported[kind])
    assert (n_layers, n_filters) == (N_DIL, FILTERS[kind])
    values = dataset_values(exported[kind])
    assert set(reads) == set(values)   # every weight in the file, and nothing missing
    for name, value in reads.items():
        layer, _, weight = name.split("/")
        want = models[kind].get_layer(layer).get_weights()[{"kernel:0": 0, "bias:0": 1}[weight]]
        np.testing.assert_array_equal(value, want)


def test_bpnetlite_needs_the_add_renumbering(models, tmp_path):
    # the no-bias model's Keras 3 Add names (add_4 ... add_7) would make bpnet-lite count 7 dilated layers
    path = export_legacy_h5(models["nobias"], tmp_path / "raw_names.h5", normalize_names=False)
    with pytest.raises(KeyError):
        bpnetlite_from_chrombpnet_reads(path)


def test_bpnetlite_reader_on_the_1x_bias_file(goldens_dir):
    # the reimplementation reads the real 1.x file the way bpnet-lite does
    path = trace_file(goldens_dir, "trace_bias_128x4")
    n_layers, n_filters, reads = bpnetlite_from_chrombpnet_reads(path)
    assert (n_layers, n_filters) == (4, 128)
    assert set(reads) == set(dataset_values(path))


# ---- (4) the 1.x reference files, re-exported ----

@pytest.mark.parametrize("which", ["init.h5", "final.h5"])
@pytest.mark.parametrize("trace", TRACES)
def test_reexported_1x_files_are_identical(trace, which, goldens_dir, tmp_path):
    """Loading a TF-Keras 2.12 file and exporting it gives the same file, except that there is no
    training_config / optimizer_weights: the same model_config text (count-head Lambda included), attributes,
    datasets."""
    golden = trace_file(goldens_dir, trace, which)
    model = model_io.load_model(golden)
    path = export_legacy_h5(model, tmp_path / "reexported.h5")

    golden_attrs, golden_keys = root_attrs(golden)
    attrs, keys = root_attrs(path)
    assert attrs == {k: v for k, v in golden_attrs.items() if k != "training_config"}
    assert keys == [k for k in golden_keys if k != "optimizer_weights"]
    with h5py.File(golden, "r") as g, h5py.File(path, "r") as f:
        assert [f.attrs[k] for k in ("keras_version", "backend")] == [g.attrs[k] for k in ("keras_version", "backend")]
    assert weights_tree(path) == weights_tree(golden)   # names, attributes (and HDF5 types), shapes, dtypes
    got, want = dataset_values(path), dataset_values(golden)
    assert list(got) == list(want)
    for name in want:
        np.testing.assert_array_equal(got[name], want[name], err_msg=name)

    assert config_diff(model_config(path), model_config(golden)) == []
    with h5py.File(golden, "r") as g, h5py.File(path, "r") as f:
        assert f.attrs["model_config"] == g.attrs["model_config"]   # the JSON text, byte for byte
    x = one_hot()
    for a, b in zip(predict(model_io.load_model(path), x), predict(model, x)):
        np.testing.assert_allclose(a, b, rtol=0, atol=1e-6)


# ---- CLI ----

def test_cli_parses():
    args = parsers.read_parser(["export", "-m", "m.h5", "-o", "out.h5"])
    assert (args.cmd, args.model_h5, args.output, args.format, args.count_head) == \
        ("export", "m.h5", "out.h5", "legacy-h5", "bytecode")
    assert parsers.read_parser(["export", "-m", "m.keras", "-o", "o.h5", "--legacy-h5"]).format == "legacy-h5"
    assert parsers.read_parser(["export", "-m", "m.h5", "-o", "o.h5", "--count-head", "named"]).count_head == "named"
    with pytest.raises(SystemExit):
        parsers.read_parser(["export", "-m", "m.h5", "-o", "o.h5", "--count-head", "source"])


def test_cli_help(capsys):
    with pytest.raises(SystemExit) as exc:
        parsers.read_parser(["export", "--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "--legacy-h5" in out and "--count-head {bytecode,named}" in out


@pytest.mark.parametrize("flags,lambda_config", [([], LAMBDA_1X), (["--count-head", "named"], LAMBDA_NAMED)])
def test_cli_export(flags, lambda_config, models, tmp_path, monkeypatch):
    import chrombpnet.CHROMBPNET as cli
    out = tmp_path / "bias.legacy.h5"
    monkeypatch.setattr(sys, "argv", ["chrombpnet", "export", "-m", models["paths"]["bias"], "-o", str(out)] + flags)
    cli.main()
    assert weights_tree(out)[""][1]["layer_names"] == BIAS_LAYER_NAMES
    out = tmp_path / "chrombpnet.legacy.h5"
    monkeypatch.setattr(sys, "argv", ["chrombpnet", "export", "-m", models["paths"]["full"], "-o", str(out)] + flags)
    cli.main()
    assert lambda_layer(model_config(out))[0]["config"] == lambda_config


# ---- TF-Keras itself (optional) ----

# argv: x.npy, then (mode, export, reference) triples. Modes: "plain" loads with no custom objects at all, "legacy"
# the way chrombpnet 1.x and the variant-scorer load models (their load_model_wrapper), "named" with the function of
# the named-function count head. The reference file holds the Keras 3 outputs (of the model and of its nested
# models) and every source weight at its TF-Keras path, <top layer>/<layer>/<weight>:0.
TF_CHECK = textwrap.dedent("""
    import sys
    import h5py
    import numpy as np
    import tensorflow as tf

    def chrombpnet_logsumexp(x):
        return tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)

    def multinomial_nll(true_counts, logits):
        raise NotImplementedError

    def load(mode, path):
        if mode == "plain":
            return tf.keras.models.load_model(path, compile=False)
        if mode == "legacy":   # chrombpnet 1.x / variant-scorer load_model_wrapper
            tf.keras.utils.get_custom_objects().update({"multinomial_nll": multinomial_nll, "tf": tf})
            return tf.keras.models.load_model(path, compile=False)
        assert mode == "named", mode
        return tf.keras.models.load_model(
            path, compile=False, custom_objects={"chrombpnet_logsumexp": chrombpnet_logsumexp})

    def check_outputs(model, x, ref, what):
        got = model.predict(x, verbose=0)
        for out, name in zip(got, ("profile", "counts")):
            np.testing.assert_allclose(out, ref[name][()], rtol=0, atol=1e-4, err_msg=what + " " + name)

    def n_datasets(group):
        found = []
        group.visititems(lambda name, obj: found.append(name) if isinstance(obj, h5py.Dataset) else None)
        return len(found)

    x = np.load(sys.argv[1])
    runs = [sys.argv[i:i + 3] for i in range(2, len(sys.argv), 3)]
    for mode in ("plain", "legacy", "named"):   # "legacy" registers custom objects globally: plain loads first
        for _, path, ref_path in [run for run in runs if run[0] == mode]:
            model = load(mode, path)
            with h5py.File(path, "r") as f, h5py.File(ref_path, "r") as ref:
                # each variable, by its path in the file, equals that dataset and the source weight of that path
                n = 0
                for top in model.layers:
                    for sub in (top.layers if isinstance(top, tf.keras.Model) else [top]):
                        for v in sub.weights:
                            key = "{}/{}/{}".format(top.name, sub.name, v.name.split("/")[-1])
                            np.testing.assert_array_equal(v.numpy(), f["model_weights"][key][()], err_msg=key)
                            np.testing.assert_array_equal(v.numpy(), ref["weights"][key][()], err_msg=key)
                            n += 1
                assert n == n_datasets(f["model_weights"]) == n_datasets(ref["weights"]), path
                for name in ref.attrs["frozen"]:
                    assert not model.get_layer(name).trainable_weights, (path, name)
                assert len(model.trainable_weights) == ref.attrs["n_trainable"], path
                check_outputs(model, x, ref["outputs"], path)
                for name in ref["nested_outputs"]:
                    check_outputs(model.get_layer(name), x, ref["nested_outputs"][name], path + " " + name)
            print("ok", mode, path)
""")


def write_tf_reference(path, kind, models, x):
    """The reference file TF_CHECK compares a TF-Keras-loaded export of `kind` with."""
    source = models["full" if kind == "full_named" else kind]
    nested = {"model_wo_bias": models["nobias"], "model": models["bias"]} if source is models["full"] else {}
    with h5py.File(path, "w") as f:
        for group, model in [("outputs", source)] + [("nested_outputs/" + n, m) for n, m in nested.items()]:
            profile, counts = predict(model, x)
            f[group + "/profile"], f[group + "/counts"] = profile, counts
        f.require_group("nested_outputs")
        for key, value in source_weights("full" if nested else kind, models).items():
            f["weights/" + key] = value
        f.attrs["frozen"] = ["model"] if nested else []
        f.attrs["n_trainable"] = len(source.trainable_weights)


def test_tf_keras_loads_the_exports(models, exported, tmp_path):
    python = os.environ.get("CHROMBPNET_TF_PYTHON")
    if not python:
        pytest.skip("CHROMBPNET_TF_PYTHON (a Python with TensorFlow 2.x) is not set")
    x = one_hot()
    np.save(tmp_path / "x.npy", x)
    argv = [python, "-c", TF_CHECK, str(tmp_path / "x.npy")]
    # the default exports load with no custom objects, and the full model as the legacy readers load it
    runs = [("plain", "bias"), ("plain", "nobias"), ("plain", "full"), ("legacy", "full"), ("named", "full_named")]
    for mode, kind in runs:
        reference = str(tmp_path / (kind + ".ref.h5"))
        if not os.path.exists(reference):
            write_tf_reference(reference, kind, models, x)
        argv += [mode, exported[kind], reference]
    result = subprocess.run(argv, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr[-3000:]
    assert result.stdout.count("ok ") == len(runs), result.stdout
