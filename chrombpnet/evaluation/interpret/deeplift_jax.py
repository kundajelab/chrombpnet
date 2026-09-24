"""DeepLIFT / DeepSHAP rules in JAX, reproducing chrombpnet 1.x (kundajelab-shap ``TFDeepExplainer``).

Batch layout: a joint batch ``[x_1 .. x_P, r_1 .. r_P]`` of 2P rows in which row p (an input) pairs with row
P + p (its reference). Every model op is row-wise, so pairs never interact.

Legacy semantics (kundajelab-shap==1, shap/explainers/deep/deep_tf.py):
* ReLU: rescale rule, ``m = (relu(zx) - relu(zr)) / (zx - zr)``, falling back to the plain gradient
  ``(z > 0)`` where ``|zx - zr| < 1e-6`` (nonlinearity_1d_handler, deep_tf.py:480-498). Every other op in a
  BPNet (conv, bias, crop, add, flatten, mean, dense) is linear and keeps its ordinary gradient.
* Counts target: ``sum(outputs[1], -1)`` (chrombpnet interpret.py).
* Profile target: ``sum(softmax(stop_gradient(mn)) * mn)`` with mean-normalised logits ``mn``. deep_tf.py:595
  (``op_handlers["StopGradient"] = break_dependence``) is commented out, so the Mul sees two varying inputs and
  gets the two-input Shapley rule (nonlinearity_2d_handler, deep_tf.py:500-528). Its logits branch is
  ``0.5 * (s_x + s_r)`` (zeroed where ``|mn_x - mn_r| < 1e-7``); the softmax branch dies at the StopGradient. So
  the logit weights are the MIDPOINT of the input's and the reference's softmax, not ``softmax(mn_x)``.
* Only the input half's multipliers are used; hypothetical scores are
  ``mean_k(m_k - sum_c r_k * m_k)`` (shap_utils.combine_mult_and_diffref).
"""
import jax
import jax.numpy as jnp

# deep_tf.py:494 (nonlinearity_1d_handler): |delta_in| < 1e-6 (strict) -> plain gradient
EPS_RELU = 1e-6
# deep_tf.py:525-526 (nonlinearity_2d_handler): |delta_in| < 1e-7 (strict) -> zero
EPS_MUL = 1e-7

PROFILE_WEIGHTINGS = ("chrombpnet", "chrombpnet_tf", "softmax_x", "tangermeme")


def split_pairs(a):
    """Split a joint batch into its input half and its reference half."""
    p = a.shape[0] // 2
    return a[:p], a[p:]


def rescale_multipliers(zx, zr):
    """ReLU rescale-rule multipliers for the input half and the reference half.

    The secant is shared by both halves; only the |dz| < EPS_RELU fallback (each half's own plain gradient)
    differs between them.
    """
    delta = zx - zr
    small = jnp.abs(delta) < EPS_RELU
    secant = (jnp.maximum(zx, 0) - jnp.maximum(zr, 0)) / jnp.where(small, jnp.ones_like(delta), delta)
    mx = jnp.where(small, (zx > 0).astype(zx.dtype), secant)
    mr = jnp.where(small, (zr > 0).astype(zr.dtype), secant)
    return mx, mr


@jax.custom_vjp
def dl_relu(z):
    """ReLU whose gradient is the DeepLIFT rescale multiplier on the input half and zero on the reference half.

    Zeroing the reference half is exact for chrombpnet targets: they read reference rows only through
    stop_gradient, so the reference half's cotangent is zero anyway and never mixes with the input half.
    """
    return jnp.maximum(z, 0)


def _dl_relu_fwd(z):
    zx, zr = split_pairs(z)
    mx, _ = rescale_multipliers(zx, zr)
    return jnp.maximum(z, 0), mx


def _dl_relu_bwd(mx, g):
    gx, gr = split_pairs(g)
    return (jnp.concatenate([gx * mx, jnp.zeros_like(gr)], axis=0),)


dl_relu.defvjp(_dl_relu_fwd, _dl_relu_bwd)


@jax.custom_vjp
def dl_relu_full(z):
    """Both-halves variant, as TFDeepExplainer computes it (reference half gets its own fallback). For tests."""
    return jnp.maximum(z, 0)


def _dl_relu_full_fwd(z):
    zx, zr = split_pairs(z)
    mx, mr = rescale_multipliers(zx, zr)
    return jnp.maximum(z, 0), jnp.concatenate([mx, mr], axis=0)


def _dl_relu_full_bwd(m, g):
    return (g * m,)


dl_relu_full.defvjp(_dl_relu_full_fwd, _dl_relu_full_bwd)


def mean_normalize(logits):
    """Mean-normalised profile logits, (rows, outlen). "Adjustments for Softmax Layers" in the DeepLIFT paper."""
    logits = jnp.reshape(logits, (logits.shape[0], -1))
    return logits - jnp.mean(logits, axis=1, keepdims=True)


def profile_weights(mn_x, mn_r, mode="chrombpnet"):
    """Per-logit weights w (pairs, outlen) such that the profile multipliers are the gradient of sum(w * mn_x).

    chrombpnet     legacy midpoint rule 0.5 * (softmax(mn_x) + softmax(mn_r)), zero where |mn_x - mn_r| < 1e-7
    chrombpnet_tf  same rule evaluated with the literal float formula of deep_tf.py:515-523 (diagnostics)
    softmax_x      plain autodiff of stop_gradient(softmax(mn)) * mn, i.e. softmax(mn_x) (comparison only)
    tangermeme     bpnet-lite / tangermeme secant of mn * softmax(mn), (y_x - y_r) / (mn_x - mn_r); where
                   |mn_x - mn_r| < 1e-6 this uses softmax(mn_x), an approximation of torch's plain gradient
                   (comparison only)
    """
    s_x = jax.nn.softmax(mn_x, axis=-1)
    s_r = jax.nn.softmax(mn_r, axis=-1)
    delta = mn_x - mn_r
    if mode == "chrombpnet":
        return jnp.where(jnp.abs(delta) < EPS_MUL, jnp.zeros_like(s_x), 0.5 * (s_x + s_r))
    if mode == "chrombpnet_tf":
        small = jnp.abs(delta) < EPS_MUL
        # out1 = 0.5 * (out11 - out10 + out01 - out00) / delta_in1 with out_ab = x0_a * x1_b, x0 = softmax, x1 = mn
        out1 = 0.5 * (s_x * mn_x - s_x * mn_r + s_r * mn_x - s_r * mn_r)
        return jnp.where(small, jnp.zeros_like(s_x), out1 / jnp.where(small, jnp.ones_like(delta), delta))
    if mode == "softmax_x":
        return s_x
    if mode == "tangermeme":
        small = jnp.abs(delta) < EPS_RELU
        secant = (mn_x * s_x - mn_r * s_r) / jnp.where(small, jnp.ones_like(delta), delta)
        return jnp.where(small, s_x, secant)
    raise ValueError("profile_weighting must be one of {}, got {!r}".format(PROFILE_WEIGHTINGS, mode))


def targets(logits, logcounts, mode="chrombpnet", full=False):
    """Scalar DeepSHAP targets of a joint batch plus the per-pair linearised deltas they must add up to.

    Returns ((T_counts, T_profile), aux) with aux = {"counts": c(x) - c(r), "profile": sum(w * (mn_x - mn_r))}
    per pair. With full=True the targets also sum the reference rows (as TFDeepExplainer does); the input half's
    multipliers are the same either way.
    """
    counts = jnp.sum(jnp.reshape(logcounts, (logcounts.shape[0], -1)), axis=-1)
    mn = mean_normalize(logits)
    mn_x, mn_r = split_pairs(jax.lax.stop_gradient(mn))
    w = jax.lax.stop_gradient(profile_weights(mn_x, mn_r, mode))
    mn_xg, mn_rg = split_pairs(mn)
    c_x, c_r = split_pairs(counts)
    if full:
        t_counts = jnp.sum(counts)
        t_profile = jnp.sum(w * mn_xg) + jnp.sum(w * mn_rg)
    else:
        t_counts = jnp.sum(c_x)
        t_profile = jnp.sum(w * mn_xg)
    aux = {
        "counts": jax.lax.stop_gradient(c_x - c_r),
        "profile": jnp.sum(w * (mn_x - mn_r), axis=-1),
    }
    return (t_counts, t_profile), aux


def hypothetical(mult, refs):
    """Hypothetical contributions, (S, K, L, 4) multipliers and references -> (S, L, 4).

    Same as shap_utils.combine_mult_and_diffref: mean over references of m - sum_c(r * m).
    """
    return jnp.mean(mult - jnp.sum(refs * mult, axis=-1, keepdims=True), axis=1)


def summation_delta(mult, x, refs):
    """sum_{l,c} m * (x - r) per pair: (P, L, 4) arrays -> (P,)."""
    return jnp.sum(mult * (x - refs), axis=(1, 2))
