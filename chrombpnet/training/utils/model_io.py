"""Single entry point for loading ChromBPNet models (chrombpnet 1.x TF-Keras .h5 files and Keras 3 files)."""
import io
import json
import warnings

import keras

from chrombpnet.training.utils.layers import LogSumExp, LogSumExpCompat
from chrombpnet.training.utils.losses import multinomial_nll

_NESTED_MODEL_CLASSES = ("Functional", "Model", "Sequential")


def require_jax_backend():
    """Fail early with a clear message if Keras was initialised with another backend."""
    backend = keras.backend.backend()
    if backend != "jax":
        raise RuntimeError(
            "chrombpnet needs the Keras JAX backend but Keras is using {!r}. Set KERAS_BACKEND=jax before "
            "anything imports keras (importing chrombpnet first does this).".format(backend))


def custom_objects():
    # "LogSumExpCompat": files re-saved from a loaded 1.x model before LogSumExpCompat.from_config returned a
    # plain LogSumExp
    return {"Lambda": LogSumExpCompat, "LogSumExp": LogSumExp, "LogSumExpCompat": LogSumExp,
            "multinomial_nll": multinomial_nll}


def _had_construction_node(layer_config):
    """True if TF-Keras 2.x gave this nested model a construction node 0 (TF's `_should_skip_first_node`).

    Functional models and Sequential models built with an input shape start with that node; their config lists
    an InputLayer first. A Sequential built without an input shape has no such node.
    """
    if layer_config.get("class_name") not in _NESTED_MODEL_CLASSES:
        return False
    config = layer_config.get("config", {})
    layers = config.get("layers", []) if isinstance(config, dict) else config
    return bool(layers) and isinstance(layers[0], dict) and layers[0].get("class_name") == "InputLayer"


def _is_node_ref(ref):
    return (isinstance(ref, list) and len(ref) in (3, 4) and isinstance(ref[0], str)
            and isinstance(ref[1], int) and not isinstance(ref[1], bool))


def _shift_refs(refs, names):
    """Decrement the node index of every [name, node_index, tensor_index(, kwargs)] ref to one of `names`
    inside the (possibly nested) list / dict structure `refs`. Returns True if anything changed."""
    if _is_node_ref(refs):
        if refs[0] in names and refs[1] >= 1:
            refs[1] -= 1
            return True
        return False
    if isinstance(refs, dict):
        refs = list(refs.values())
    changed = False
    for value in refs if isinstance(refs, list) else ():
        changed |= _shift_refs(value, names)
    return changed


def _shift_nested_node_indices(model_config):
    """Renumber calls to nested models in a TF-Keras 2.x functional config for Keras 3.

    TF-Keras counts a nested model's own construction as its inbound node 0, so the first call inside the outer
    model is node 1 (e.g. `["model_wo_bias", 1, 0, {}]` in chrombpnet.h5). Keras 3 does not create that node.
    Keras 3 already subtracts 1 from input_layers / output_layers refs to a Functional (functional_from_config's
    get_tensor), but not from the legacy list-form inbound_nodes refs, nor from any ref to a Sequential. So:
    inbound_nodes refs are shifted for every nested model that had a construction node, input_layers /
    output_layers refs only for such Sequential models. Returns True if anything changed.
    """
    config = model_config.get("config", {})
    if not isinstance(config, dict):
        return False
    layers = config.get("layers", [])
    with_node0 = {layer["name"]: layer.get("class_name") for layer in layers if _had_construction_node(layer)}
    changed = False
    for layer in layers:
        if layer.get("class_name") in _NESTED_MODEL_CLASSES:
            changed |= _shift_nested_node_indices(layer)
        for node in layer.get("inbound_nodes", []):
            # legacy list form only: Keras 3 dict-form nodes already hold Keras 3 node indices
            if isinstance(node, list):
                changed |= _shift_refs(node, set(with_node0))
    sequential = {name for name, class_name in with_node0.items() if class_name == "Sequential"}
    for key in ("input_layers", "output_layers"):
        changed |= _shift_refs(config.get(key, []), sequential)
    return changed


def _load_legacy_h5(path, compile):
    """Load a TF-Keras 2.x .h5 whose graph calls nested models, fixing the config on an in-memory copy."""
    import h5py
    from keras.src.legacy.saving import legacy_h5_format

    with open(path, "rb") as fh:
        buffer = io.BytesIO(fh.read())
    with h5py.File(buffer, "r+") as f:
        config = json.loads(f.attrs["model_config"])
        if _shift_nested_node_indices(config):
            f.attrs["model_config"] = json.dumps(config)
        return legacy_h5_format.load_model_from_hdf5(f, custom_objects=custom_objects(), compile=compile)


def _is_legacy_h5_with_nested_models(path):
    if not str(path).endswith((".h5", ".hdf5")):
        return False
    import h5py

    with h5py.File(path, "r") as f:
        version = f.attrs.get("keras_version", "")
        version = version.decode() if isinstance(version, bytes) else str(version)
        if not version.startswith("2.") or "model_config" not in f.attrs:
            return False
        config = json.loads(f.attrs["model_config"])
    return any(layer.get("class_name") in _NESTED_MODEL_CLASSES
               for layer in config.get("config", {}).get("layers", []))


def load_model(path, compile=False):
    """Load a bias / chrombpnet / chrombpnet_nobias model.

    Works for .keras files, .h5 files written by Keras 3, and legacy TF-Keras 2.x .h5 files, including full
    `chrombpnet.h5` models (nested bias + no-bias models and a logsumexp Lambda count head). Optimizer state
    is never restored.
    """
    if keras.backend.backend() != "jax":
        warnings.warn("Loading {} with the Keras {!r} backend; chrombpnet is tested on JAX only.".format(
            path, keras.backend.backend()))
    if _is_legacy_h5_with_nested_models(path):
        return _load_legacy_h5(str(path), compile)
    return keras.saving.load_model(str(path), compile=compile, custom_objects=custom_objects())


def load_model_wrapper(model_h5):
    """Backwards-compatible name used by chrombpnet 1.x code and downstream pipelines."""
    model = load_model(model_h5)
    print("got the model")
    model.summary()
    return model
