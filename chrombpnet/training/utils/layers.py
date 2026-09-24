import base64

import keras
from keras import ops


@keras.saving.register_keras_serializable(package="chrombpnet")
class LogSumExp(keras.layers.Layer):
    """Combine the bias and no-bias log-count heads: logsumexp over the last axis, keepdims.

    Replaces the `Lambda(lambda x: tf.math.reduce_logsumexp(...))` of chrombpnet 1.x, which cannot be
    serialized portably (marshalled Python bytecode) and cannot run on the JAX backend.
    """

    def call(self, inputs):
        return ops.logsumexp(inputs, axis=-1, keepdims=True)

    def compute_output_shape(self, input_shape):
        return tuple(input_shape[:-1]) + (1,)


_LAMBDA_ONLY_KEYS = ("function", "output_shape", "mask", "arguments", "module", "function_type",
                     "output_shape_type", "output_shape_module")

# Name of the function the count-head Lambda of `chrombpnet export --legacy-h5 --count-head named` files refers to
# (function_type "function"): TF-Keras 2.x resolves it through custom_objects, chrombpnet maps it to LogSumExp.
LOGSUMEXP_LAMBDA_FUNCTION = "chrombpnet_logsumexp"


def _lambda_code_bytes(function):
    """Raw marshalled code of a serialized Lambda function (list or dict form), or b"" if unknown."""
    code = None
    if isinstance(function, (list, tuple)) and function:
        code = function[0]
    elif isinstance(function, dict):
        code = function.get("config", {}).get("code", function.get("code"))
    elif isinstance(function, str):
        code = function
    if not isinstance(code, str):
        return b""
    try:
        return base64.decodebytes(code.encode("ascii"))
    except (ValueError, TypeError):
        return b""


def _is_logsumexp_lambda(config):
    """The 1.x bytecode Lambda (`tf.math.reduce_logsumexp`) or the named-function Lambda of legacy exports.

    Keras 3's legacy loader drops `function_type` / `module` before from_config, so a named function arrives as
    just the name.
    """
    function = config.get("function")
    if config.get("function_type", "function") == "function" and function == LOGSUMEXP_LAMBDA_FUNCTION:
        return True
    return b"reduce_logsumexp" in _lambda_code_bytes(function)


class LogSumExpCompat(LogSumExp):
    """Stand-in for the logsumexp Lambda inside chrombpnet 1.x `chrombpnet.h5` files.

    Passed to load_model as custom_objects={"Lambda": LogSumExpCompat}. It never unmarshals the stored bytecode
    (which only loads on the Python version that wrote it). `chrombpnet export --legacy-h5` files carry the same
    Lambda, or with `--count-head named` one that names the function (function_type "function", function
    LOGSUMEXP_LAMBDA_FUNCTION) instead of storing bytecode, which is accepted too. Any other Lambda is refused
    rather than silently replaced. from_config returns a plain LogSumExp, so a loaded 1.x model re-saves as
    `chrombpnet>LogSumExp` and loads again.
    """

    @classmethod
    def from_config(cls, config):
        config = dict(config)
        name = config.get("name")
        if name != "logcount_predictions" or not _is_logsumexp_lambda(config):
            raise ValueError(
                "Cannot load Lambda layer {!r}: only the logsumexp count head of chrombpnet 1.x models "
                "(`logcount_predictions`) is supported. Rebuild custom architectures without Lambda layers "
                "(see chrombpnet.training.utils.layers).".format(name))
        for key in _LAMBDA_ONLY_KEYS:
            config.pop(key, None)
        return LogSumExp(**config)
