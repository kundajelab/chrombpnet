"""multinomial_nll against a float64 reference of tfp's Multinomial(total_count=sum(y), logits).log_prob(y)."""
import numpy as np
import pytest
from scipy.special import gammaln, log_softmax

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)

from chrombpnet.training.utils.losses import multinomial_nll

BINS = 1000


def reference_nll(y, logits):
    """-(sum_i y_i log p_i + lgamma(n + 1) - sum_i lgamma(y_i + 1)), n = sum_i y_i, in float64."""
    y, logits = np.asarray(y, np.float64), np.asarray(logits, np.float64)
    weighted = np.where(y > 0, y * log_softmax(logits, axis=-1), 0.0)
    return -(weighted.sum(-1) + gammaln(y.sum(-1) + 1) - gammaln(y + 1).sum(-1))


def check(y, logits):
    got = np.asarray(multinomial_nll(y.astype(np.float32), logits.astype(np.float32)))
    want = reference_nll(y, logits)
    assert got.shape == want.shape == (y.shape[0],)
    # float32 over 1000 bins: ~1e-3 absolute on losses of a few hundred
    np.testing.assert_allclose(got, want, rtol=1e-5, atol=5e-3)


@pytest.fixture
def logits():
    return np.random.RandomState(0).normal(size=(6, BINS))


def test_integer_counts(logits):
    rng = np.random.RandomState(1)
    y = rng.poisson(0.5, (6, BINS)).astype(np.float64)
    y[0] = 0          # an empty profile (n = 0)
    y[1] = 0
    y[1, 7] = 1       # n = 1
    check(y, logits)


def test_fractional_counts(logits):
    # e.g. a normalized bigwig: most bins are fractions in (0, 1), for which lgamma(y + 1) != 0
    rng = np.random.RandomState(2)
    y = rng.poisson(0.5, (6, BINS)) * 0.37
    y[0] = 0
    y[0, :2] = 0.25   # total n in (0, 1)
    y[1] = 0
    y[1, :4] = 0.25   # total n = 1 exactly, from fractional bins
    check(y, logits)


def test_zero_counts_with_minus_inf_logits():
    # a zero count contributes 0 even where log_softmax is -inf
    y = np.array([[0.0, 3.0, 1.0]])
    logits = np.array([[-np.inf, 0.0, 1.0]])
    got = np.asarray(multinomial_nll(y.astype(np.float32), logits.astype(np.float32)))
    assert np.isfinite(got).all()
    np.testing.assert_allclose(got, reference_nll(y, logits), rtol=1e-6)
