import keras
from keras import ops


def _gammaln(x):
    """log|Gamma(x)|; keras.ops has no lgamma, so dispatch on the active backend."""
    backend = keras.backend.backend()
    if backend == "jax":
        from jax.scipy.special import gammaln
        return gammaln(x)
    if backend == "tensorflow":
        import tensorflow as tf
        return tf.math.lgamma(x)
    if backend == "torch":
        import torch
        return torch.lgamma(x)
    from scipy.special import gammaln
    return gammaln(x)


#from https://github.com/kundajelab/basepair/blob/cda0875571066343cdf90aed031f7c51714d991a/basepair/losses.py#L87
@keras.saving.register_keras_serializable(package="chrombpnet")
def multinomial_nll(true_counts, logits):
    """Compute the multinomial negative log-likelihood
    Args:
      true_counts: observed count values
      logits: predicted logit values

    Same value as the tensorflow_probability Multinomial(total_count, logits).log_prob used by chrombpnet 1.x
    (including the count-only log-factorial terms, so logged losses stay comparable). Returns one value per
    example; Keras averages over the batch, which equals the old pre-averaged scalar.
    """
    true_counts = ops.cast(true_counts, "float32")
    logits = ops.cast(logits, "float32")
    counts_per_example = ops.sum(true_counts, axis=-1)
    log_probs = ops.log_softmax(logits, axis=-1)
    # multiply_no_nan: a zero count contributes 0 even where log_prob is -inf
    weighted = ops.where(true_counts > 0, true_counts * log_probs, 0.0)
    # log(n!) and sum_i log(y_i!) as lgamma(c + 1); lgamma(1) = lgamma(2) = 0 exactly, so skip c = 0 and c = 1
    # (float32 lgamma is not exact there and ~1000 mostly-zero bins would otherwise add a spurious offset).
    # Fractional counts (e.g. a normalized bigwig) keep their nonzero lgamma(c + 1), as in tfp.
    log_factorial = lambda c: ops.where(ops.logical_or(ops.equal(c, 0.0), ops.equal(c, 1.0)), 0.0,
                                        _gammaln(c + 1.0))
    log_likelihood = (ops.sum(weighted, axis=-1)
                      + log_factorial(counts_per_example)
                      - ops.sum(log_factorial(true_counts), axis=-1))
    return -log_likelihood
