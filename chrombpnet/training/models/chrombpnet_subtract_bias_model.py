import tensorflow as tf
from tensorflow.keras.layers import Input, Subtract
from tensorflow.keras.models import Model, load_model
from tensorflow.keras.utils import get_custom_objects
from chrombpnet.training.utils.losses import multinomial_nll


def _load_frozen(path, name):
    get_custom_objects().update({"multinomial_nll": multinomial_nll, "tf": tf})
    m = load_model(path, compile=False)
    m._name = name
    for L in m.layers:
        L.trainable = False
    return m


def build_subtract_model(bpnet_h5, scaled_bias_h5, output_h5):
    """
    Build a Keras wrapper that subtracts a (frozen, scaled) bias model from a
    (frozen) BPNet model and saves it as `output_h5`.

    Output semantics:
      - profile head: `bpnet_logits - bias_logits`. The downstream multinomial
        NLL is invariant to additive constants per position, so subtracting a
        bias profile gives the bias-corrected profile-shape logits.
      - counts head: `bpnet_log_counts - bias_log_counts` (a log-count *delta*,
        not log of the count delta). Both heads are now linear in the two
        sub-models' outputs, so DeepSHAP / TFDeepExplainer traces through cleanly:
            attr_nobias(x) = attr_bpnet(x) - attr_bias(x)
        for both profile and counts.

    The earlier `log(max(exp(bp) - exp(bias), eps))` formulation produced a
    bias-corrected count in the count space, but introduced a non-linear Lambda
    (with a hard floor) that broke gradient flow for SHAP whenever the bias
    matched or exceeded the BPNet prediction. The log-count-delta form keeps
    the Lambda out and the gradients linear at the cost of changing the count
    interpretation: downstream code that expects `exp(count_out)` to be a real
    cut-count must translate (the delta is between two log-counts, not a log
    of a count).
    """
    bpnet = _load_frozen(bpnet_h5, "bpnet")
    bias = _load_frozen(scaled_bias_h5, "scaled_bias")

    sequence_len = int(bpnet.input_shape[1])
    inp = Input(shape=(sequence_len, 4), name='sequence')

    bp_p, bp_lc = bpnet(inp)
    bs_p, bs_lc = bias(inp)

    profile_out = Subtract(name='logits_profile_predictions')([bp_p, bs_p])
    count_out = Subtract(name='logcount_predictions')([bp_lc, bs_lc])

    model = Model(inputs=[inp], outputs=[profile_out, count_out])
    model.save(output_h5)
    return model
