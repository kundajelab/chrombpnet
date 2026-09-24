"""Export a ChromBPNet model as a TF-Keras 2.x full-model .h5 file (`chrombpnet export --legacy-h5`).

Keras 3 writes .h5 files that TF-Keras 2.x (chrombpnet 1.x, the kundajelab variant-scorer) cannot load, and whose
weight datasets are named `layer/kernel` instead of `layer/kernel:0` (so h5py readers such as bpnet-lite's
`BPNet.from_chrombpnet` / `ChromBPNet.from_chrombpnet` do not find them). This module writes the layout TF-Keras
2.12 wrote for chrombpnet 1.x models:

* root attributes `keras_version` ("2.12.0"), `backend` ("tensorflow") and `model_config`: the Keras 2 functional
  JSON (class_name "Functional", list-form inbound_nodes, InputLayer `batch_input_shape`, plain dtype strings,
  no `module` / `registered_name` / `build_config` keys). A nested model called once is referenced with node
  index 1, as TF-Keras counts its construction as node 0.
* `model_weights` with `layer_names`, `backend` and `keras_version` attributes, one group per layer holding a
  `weight_names` attribute (`bpnet_1conv/kernel:0`, ...) and the datasets
  `model_weights/<layer>/<layer>/kernel:0`; a nested model's group holds all its weights,
  `model_weights/<nested model>/<layer>/kernel:0`, in TF's order (the trainable weights of its layers, then the
  non-trainable ones), and an empty `top_level_model_weights` group.

There is no `training_config` and no `optimizer_weights`: the file is meant for `load_model(..., compile=False)`.

Auto-generated Keras 3 names are replaced by the ones TF-Keras gave the same layers in chrombpnet 1.x: a
Functional named `functional` / `functional_N` becomes `model` (`model_1`, ... for a second one in the same
model) and weightless layers named after their class (`add_4`, `add_5`, ...) are renumbered from scratch per
model (`add`, `add_1`, ...), unless that would clash with another layer name (then all names are kept, with a
warning). Layers with weights keep their names, so every weight keeps its name. The renumbering matters to
bpnet-lite, which takes the number of dilated layers from the largest number at the end of any layer name
(`add_7` in a 4-layer model would make it look for `bpnet_7conv`).

The count head of a full chrombpnet model (the registered `LogSumExp` layer) is a Keras 2 `Lambda`, in one of two
forms (`count_head`, CLI `--count-head`):

* "bytecode" (the default): the `Lambda` of the chrombpnet 1.x `chrombpnet.h5` files, verbatim. It stores
  `lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)` as Python 3.8 marshalled bytecode
  (`function_type: "lambda"`, module `chrombpnet.training.models.chrombpnet_with_bias_model`). TF-Keras
  unmarshals it on Python 3.8 - 3.10 (tested: TF 2.8 on 3.9, TF 2.12 on 3.10), which covers chrombpnet 1.x and
  the kundajelab variant-scorer: their `load_model(path, compile=False)` calls (custom_objects `multinomial_nll`
  and `tf` only, or none) load the file unchanged. Python >= 3.11 cannot unmarshal it
  (`ValueError: bad marshal data`).
* "named": for TF-Keras readers on Python >= 3.11. A `Lambda` that names its function
  (`function_type: "function"`, `function: "chrombpnet_logsumexp"`), which TF-Keras resolves through
  `custom_objects` (without it, loading fails with `AttributeError: 'NoneType' object has no attribute 'get'`)::

    import tensorflow as tf

    def chrombpnet_logsumexp(x):
        return tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)

    model = tf.keras.models.load_model("chrombpnet.legacy.h5", compile=False,
                                       custom_objects={"chrombpnet_logsumexp": chrombpnet_logsumexp})

  (or `tf.keras.utils.get_custom_objects()["chrombpnet_logsumexp"] = chrombpnet_logsumexp` before loading).

Bias and no-bias models have no Lambda and load without custom objects. chrombpnet's own `model_io.load_model`
reads all of these files (both count-head forms, without unmarshalling anything). bpnet-lite reads only the weight
datasets of bias / no-bias files.
"""
import json
import os
import re
import warnings

import numpy as np

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)
import h5py
import keras
from keras.src.legacy.saving.legacy_h5_format import save_attributes_to_hdf5_group
from keras.src.models.functional import Functional
from keras.src.ops.function import make_node_key
from keras.src.utils.naming import to_snake_case

from chrombpnet.training.utils.layers import LOGSUMEXP_LAMBDA_FUNCTION, LogSumExp

LEGACY_KERAS_VERSION = "2.12.0"   # the TF-Keras version that wrote the chrombpnet 1.x reference files
LEGACY_BACKEND = "tensorflow"

# The count-head Lambda of the chrombpnet 1.x chrombpnet.h5 files (TF-Keras 2.12 on Python 3.8), verbatim:
# `lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)` as base64 marshalled Python 3.8 bytecode
LEGACY_LOGSUMEXP_BYTECODE = (
    "4wEAAAAAAAAAAAAAAAEAAAAFAAAAUwAAAHMSAAAAdABqAWoCfABkAWQCZAONA1MAKQRO6f////9U\n"
    "KQLaBGF4aXPaCGtlZXBkaW1zKQPaAnRm2gRtYXRo2hByZWR1Y2VfbG9nc3VtZXhwKQHaAXipAHII\n"
    "AAAA+kgvb3B0L2Nocm9tYnBuZXQvY2hyb21icG5ldC90cmFpbmluZy9tb2RlbHMvY2hyb21icG5l\n"
    "dF93aXRoX2JpYXNfbW9kZWwucHnaCDxsYW1iZGE+gAAAAPMAAAAA\n")
LEGACY_LOGSUMEXP_MODULE = "chrombpnet.training.models.chrombpnet_with_bias_model"
COUNT_HEADS = ("bytecode", "named")   # the default first

# Keras 2 initializers and the config keys TF-Keras 2.x accepts for them
_INITIALIZER_KEYS = {
    "GlorotUniform": ("seed",), "GlorotNormal": ("seed",), "HeUniform": ("seed",), "HeNormal": ("seed",),
    "LecunUniform": ("seed",), "LecunNormal": ("seed",), "Zeros": (), "Ones": (), "Constant": ("value",),
    "RandomNormal": ("mean", "stddev", "seed"), "RandomUniform": ("minval", "maxval", "seed"),
    "TruncatedNormal": ("mean", "stddev", "seed"), "VarianceScaling": ("scale", "mode", "distribution", "seed"),
    "Orthogonal": ("gain", "seed"), "Identity": ("gain",),
}
_CONV_KEYS = ("filters", "kernel_size", "strides", "padding", "data_format", "dilation_rate", "groups",
              "activation", "use_bias", "kernel_initializer", "bias_initializer", "kernel_regularizer",
              "bias_regularizer", "activity_regularizer", "kernel_constraint", "bias_constraint")
_DENSE_KEYS = ("units", "activation", "use_bias", "kernel_initializer", "bias_initializer", "kernel_regularizer",
               "bias_regularizer", "activity_regularizer", "kernel_constraint", "bias_constraint")
# Keras 3 class -> Keras 2 config keys after name / trainable / dtype, in the order TF-Keras 2.12 writes them
_LAYER_KEYS = {
    keras.layers.Conv1D: _CONV_KEYS,
    keras.layers.Dense: _DENSE_KEYS,
    keras.layers.Cropping1D: ("cropping",),
    keras.layers.Add: (),
    keras.layers.Concatenate: ("axis",),
    keras.layers.GlobalAveragePooling1D: ("data_format", "keepdims"),
    keras.layers.Flatten: ("data_format",),
}
# weightless layers whose Keras auto-names (class name in snake case) are renumbered per model
_RENUMBERED = (keras.layers.Add, keras.layers.Concatenate, keras.layers.Cropping1D, keras.layers.Flatten,
               keras.layers.GlobalAveragePooling1D)


def logsumexp_lambda_fields(count_head="bytecode"):
    """The Lambda config keys of the exported count head after name / trainable / dtype, in TF-Keras 2.12's
    order: the 1.x bytecode Lambda ("bytecode") or the named-function one ("named"), see the module docstring."""
    if count_head == "bytecode":
        function, function_type, module = [LEGACY_LOGSUMEXP_BYTECODE, None, None], "lambda", LEGACY_LOGSUMEXP_MODULE
    elif count_head == "named":
        function, function_type, module = LOGSUMEXP_LAMBDA_FUNCTION, "function", None
    else:
        raise ValueError("count_head must be one of {}, got {!r}.".format(", ".join(COUNT_HEADS), count_head))
    return {"function": function, "function_type": function_type, "module": module, "output_shape": None,
            "output_shape_type": "raw", "output_shape_module": None, "arguments": {}}


def _native(value):
    """JSON-ready copy: tuples -> lists, numpy scalars -> Python scalars."""
    if isinstance(value, dict):
        return {k: _native(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_native(v) for v in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def _keras2_object(serialized, what, layer):
    """A serialized Keras 3 initializer / regularizer / constraint as Keras 2 {"class_name", "config"}."""
    if serialized is None:
        return None
    if not isinstance(serialized, dict) or "class_name" not in serialized:
        raise ValueError("Layer {!r}: cannot export {} {!r} to TF-Keras 2.x.".format(layer.name, what, serialized))
    class_name, config = serialized["class_name"], dict(serialized.get("config") or {})
    if what.endswith("initializer"):
        if class_name not in _INITIALIZER_KEYS:
            raise ValueError("Layer {!r}: initializer {!r} is not supported by the legacy .h5 export.".format(
                layer.name, class_name))
        config = {k: config[k] for k in _INITIALIZER_KEYS[class_name] if k in config}
    return {"class_name": class_name, "config": _native(config)}


def _layer_dtype(layer):
    # the variable dtype: a mixed_bfloat16 layer (float32 variables) is written as a float32 layer; bfloat16
    # variables (which HDF5 / TF-Keras 2.x readers cannot take) are written as float32
    dtype = layer.dtype_policy.variable_dtype
    return "float32" if dtype == "bfloat16" else dtype


def _weight_value(variable):
    value = np.asarray(keras.ops.convert_to_numpy(variable))
    return value.astype("float32") if value.dtype.name == "bfloat16" else value


def _auto_name_base(layer):
    return to_snake_case(layer.__class__.__name__)


def _is_auto_name(name, base):
    return re.fullmatch(re.escape(base) + r"(_\d+)?", name) is not None


def _legacy_names(layers, normalize):
    """{id(layer): exported name} for the layers of one model (see the module docstring)."""
    names = {id(layer): layer.name for layer in layers}
    if not normalize:
        return names
    counters = {}
    renamed = {}
    for layer in layers:
        if isinstance(layer, Functional) and _is_auto_name(layer.name, "functional"):
            base = "model"
        elif isinstance(layer, _RENUMBERED) and not layer.weights and _is_auto_name(layer.name,
                                                                                     _auto_name_base(layer)):
            base = _auto_name_base(layer)
        else:
            continue
        n = counters.get(base, 0)
        counters[base] = n + 1
        renamed[id(layer)] = base if n == 0 else "{}_{}".format(base, n)
    candidate = {**names, **renamed}
    if len(set(candidate.values())) != len(layers):   # a new name clashes with a kept one: keep the originals
        warnings.warn("Keeping the Keras 3 layer names ({}): the TF-Keras names would clash with other layer "
                      "names. bpnet-lite counts the dilated layers from the largest number at the end of a layer "
                      "name, add_N included.".format(", ".join(
                          layer.name for layer in layers if id(layer) in renamed)))
        return names
    return candidate


def export_model_name(model, normalize=True):
    """Name the exported top-level model gets (`functional*` -> `model`)."""
    if normalize and _is_auto_name(model.name, "functional"):
        return "model"
    return model.name


class _Exporter:
    def __init__(self, normalize_names=True, count_head="bytecode"):
        self.normalize = normalize_names
        self.count_head_fields = logsumexp_lambda_fields(count_head)

    # ---- model_config ----

    def functional_config(self, model, name):
        """Keras 2 config dict of a Functional model (the value of its "config" key)."""
        if not isinstance(model, Functional):
            raise ValueError("Only functional models can be exported to a legacy .h5 file; {!r} is a {}.".format(
                model.name, type(model).__name__))
        names = _legacy_names(model.layers, self.normalize)
        layer_configs = []
        for layer in model.layers:
            nodes = []
            for index, node in enumerate(layer._inbound_nodes):
                if make_node_key(layer, index) in model._nodes and not node.is_input:
                    nodes.append(self._node_data(model, names, layer, node))
            class_name, config = self.layer_config(layer, names[id(layer)])
            layer_configs.append({"class_name": class_name, "config": config, "name": names[id(layer)],
                                  "inbound_nodes": nodes})
        return {"name": name, "trainable": bool(model.trainable), "layers": layer_configs,
                "input_layers": [self._tensor_ref(model, names, t)[:3] for t in model.inputs],
                "output_layers": [self._tensor_ref(model, names, t)[:3] for t in model.outputs]}

    @staticmethod
    def _node_index(model, operation, node_index):
        """TF-Keras node index of `operation`'s node `node_index` within `model`: only nodes of this model count,
        and a nested functional model starts at 1 (TF-Keras gave it a construction node 0)."""
        kept = 1 if isinstance(operation, Functional) else 0
        return kept + sum(1 for i in range(node_index) if make_node_key(operation, i) in model._nodes)

    def _tensor_ref(self, model, names, tensor, kwargs=None):
        operation, node_index, tensor_index = tensor._keras_history
        return [names[id(operation)], self._node_index(model, operation, node_index), tensor_index,
                {} if kwargs is None else kwargs]

    def _node_data(self, model, names, layer, node):
        args = node.arguments.args
        # Keras 3 records call arguments it filled in itself (mask=None); TF-Keras stored only explicit ones
        kwargs = {k: v for k, v in node.arguments.kwargs.items() if v is not None}
        if len(args) != 1 or any(isinstance(v, keras.KerasTensor) for v in keras.tree.flatten(kwargs)):
            raise ValueError("Layer {!r}: only calls on one tensor (or one list of tensors) can be exported.".format(
                layer.name))
        first = args[0]
        tensors = list(first) if isinstance(first, (list, tuple)) else [first]
        if not all(isinstance(t, keras.KerasTensor) for t in tensors):
            raise ValueError("Layer {!r}: only calls on tensors can be exported.".format(layer.name))
        return [self._tensor_ref(model, names, t, _native(dict(kwargs))) for t in tensors]

    def layer_config(self, layer, name):
        """(Keras 2 class_name, config) of one layer."""
        if isinstance(layer, keras.layers.InputLayer):
            config = layer.get_config()
            if config.get("optional"):
                raise ValueError("Optional inputs ({!r}) cannot be exported.".format(layer.name))
            return "InputLayer", {"batch_input_shape": _native(config["batch_shape"]), "dtype": config["dtype"],
                                  "sparse": bool(config.get("sparse", False)),
                                  "ragged": bool(config.get("ragged", False)), "name": name}
        if isinstance(layer, Functional):
            return "Functional", self.functional_config(layer, name)
        base = {"name": name, "trainable": bool(layer.trainable), "dtype": _layer_dtype(layer)}
        if isinstance(layer, LogSumExp):
            # the 1.x bytecode Lambda or a named-function one, see the module docstring
            return "Lambda", dict(base, **self.count_head_fields)
        for cls, keys in _LAYER_KEYS.items():
            if type(layer) is cls:
                break
        else:
            raise ValueError(
                "Layer {!r} ({}) is not supported by the legacy .h5 export, which covers the chrombpnet "
                "architectures: {}.".format(layer.name, type(layer).__name__, ", ".join(
                    ["InputLayer", "Functional", "LogSumExp"] + [c.__name__ for c in _LAYER_KEYS])))
        config = layer.get_config()
        if config.get("quantization_config") is not None or config.get("lora_rank"):
            raise ValueError("Layer {!r}: quantized / LoRA layers cannot be exported.".format(layer.name))
        if "activity_regularizer" in keys:
            config.setdefault("activity_regularizer", keras.regularizers.serialize(layer.activity_regularizer))
        out = dict(base)
        for key in keys:
            value = config.get(key)
            if key.endswith(("_initializer", "_regularizer", "_constraint")):
                value = _keras2_object(value, key, layer)
            elif key == "activation" and not isinstance(value, str):
                raise ValueError("Layer {!r}: only built-in activations (by name) can be exported, got {!r}.".format(
                    layer.name, value))
            out[key] = _native(value)
        return type(layer).__name__, out

    # ---- model_weights ----

    @staticmethod
    def legacy_weights(layer, trainable=True):
        """([(owner layer, variable)] trainable, [...] non-trainable) of `layer`, in TF-Keras 2.x order: a
        model lists the trainable weights of all its layers, then the non-trainable ones."""
        trainable = trainable and layer.trainable
        if isinstance(layer, Functional):
            train, frozen = [], []
            for sub in layer.layers:
                t, f = _Exporter.legacy_weights(sub, trainable)
                train += t
                frozen += f
            own = layer._trainable_variables + layer._non_trainable_variables
            if own:
                raise ValueError("Model {!r} has weights of its own, which the legacy .h5 export does not "
                                 "support.".format(layer.name))
            return train, frozen
        weights = list(layer.weights)
        if trainable:
            return ([(layer, v) for v in weights if v.trainable], [(layer, v) for v in weights if not v.trainable])
        return [], [(layer, v) for v in weights]

    def weight_entries(self, layer):
        """[(TF weight name, numpy value)] of one top-level layer. Layers with weights are never renamed, so a
        weight is named after the layer that owns it, as TF-Keras named it: `bpnet_1conv/kernel:0`."""
        train, frozen = self.legacy_weights(layer)
        return [("{}/{}:0".format(owner.name, variable.name), _weight_value(variable))
                for owner, variable in train + frozen]


def write_legacy_h5(model, out_path, normalize_names=True, count_head="bytecode"):
    """Write `model` (a Keras 3 functional model) to `out_path` in the TF-Keras 2.x full-model layout."""
    exporter = _Exporter(normalize_names, count_head)
    model_name = export_model_name(model, normalize_names)
    model_config = {"class_name": "Functional", "config": exporter.functional_config(model, model_name)}
    names = _legacy_names(model.layers, normalize_names)

    with h5py.File(out_path, "w") as f:
        # what TF-Keras 2.x save_model_to_hdf5 writes, with the same HDF5 string types: str attributes are
        # variable-length UTF-8, bytes attributes variable-length ASCII
        f.attrs["keras_version"] = LEGACY_KERAS_VERSION
        f.attrs["backend"] = LEGACY_BACKEND
        f.attrs["model_config"] = json.dumps(model_config).encode("utf8")

        weights = f.create_group("model_weights")
        save_attributes_to_hdf5_group(weights, "layer_names", [names[id(layer)].encode("utf8")
                                                               for layer in model.layers])
        weights.attrs["backend"] = LEGACY_BACKEND.encode("utf8")
        weights.attrs["keras_version"] = LEGACY_KERAS_VERSION.encode("utf8")
        # TF-Keras creates the groups sorted by layer name
        for layer in sorted(model.layers, key=lambda layer: names[id(layer)]):
            group = weights.create_group(names[id(layer)])
            entries = exporter.weight_entries(layer)
            save_attributes_to_hdf5_group(group, "weight_names", [name.encode("utf8") for name, _ in entries])
            for name, value in entries:
                dataset = group.create_dataset(name, value.shape, dtype=value.dtype)
                if value.shape:
                    dataset[:] = value
                else:
                    dataset[()] = value
        save_attributes_to_hdf5_group(weights.create_group("top_level_model_weights"), "weight_names", [])
    return out_path


def export_legacy_h5(model_or_path, out_path, normalize_names=True, count_head="bytecode"):
    """Export a bias / chrombpnet / chrombpnet_nobias model (a Keras 3 model, or a path to any .h5 / .keras file
    `model_io.load_model` reads) as a TF-Keras 2.x full-model .h5 file. The file is written next to `out_path`
    and moved into place when complete. Returns `out_path`.

    count_head: how the count head of a full chrombpnet model is written, "bytecode" (the default: the 1.x
    Lambda, which TF-Keras loads on Python 3.8 - 3.10 without custom objects) or "named" (for TF-Keras on
    Python >= 3.11, which then needs custom_objects={"chrombpnet_logsumexp": ...}); see the module docstring.
    """
    logsumexp_lambda_fields(count_head)   # refuse a bad value before loading anything
    if isinstance(model_or_path, (str, os.PathLike)):
        from chrombpnet.training.utils.model_io import load_model
        model = load_model(str(model_or_path), compile=False)
    else:
        model = model_or_path
    out_path = str(out_path)
    tmp_path = out_path + ".tmp"
    try:
        write_legacy_h5(model, tmp_path, normalize_names, count_head)
        os.replace(tmp_path, out_path)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)
    return out_path
