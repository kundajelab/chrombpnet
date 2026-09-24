"""Vendored TF-MoDISco 0.5.16 viz helpers: ic_scale numerics and logo plotting defaults.

The recorded numbers were computed with modisco 0.5.16.0's own `util.compute_per_position_ic` source
(executed in isolation) on the packaged reference motifs, after the same row renormalisation that
auto_shift_detect.ic_scale applies. They pin the Tn5/DNase shift detection inputs.
"""
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pytest  # noqa: E402

from chrombpnet.data import DefaultDataFile, get_default_data_path  # noqa: E402
from chrombpnet.helpers.preprocessing import auto_shift_detect  # noqa: E402
from chrombpnet.utils import viz_sequence  # noqa: E402

# (motif file, motif name) -> (sum of ic_scale(pwm), per-position IC of row 0), modisco 0.5.16.0
RECORDED = {
    ("atac_ref_motifs", "GSE101074_naive_hESC_ATAC_plus"): (0.9355284058615381, 0.013112558776653127),
    ("atac_ref_motifs", "scATAC_fibroblast_minus"): (0.9337631655851852, 0.011157012521581522),
    ("atac_ref_motifs", "GSE101905_d3_SSEAp_minus"): (1.0062964854733227, 0.009016269127217091),
    ("dnase_ref_motifs", "ENCSR025EYJ_DNASE_minus"): (0.8523481482696651, 0.03248084419082359),
    ("dnase_ref_motifs", "ENCSR618HFT_DNASE_plus"): (0.3119302138051853, 0.0018214787336600646),
    ("dnase_ref_motifs", "ENCSR014FPY_DNASE_plus"): (0.5520155221819234, 0.007262193885972379),
}


def ref_pwms(key):
    plus, minus = auto_shift_detect.get_ref_pwms(get_default_data_path(getattr(DefaultDataFile, key)))
    return {**plus, **minus}


def modisco_0516_ic(ppm, background, pseudocount=0.001):
    # literal transcription of modisco 0.5.16.0 util.compute_per_position_ic (without the print)
    if not np.allclose(np.sum(ppm, axis=1), 1.0, atol=1.0e-5):
        ppm = ppm / np.sum(ppm, axis=1)[:, None]
    alphabet_len = len(background)
    ic = ((np.log((ppm + pseudocount) / (1 + pseudocount * alphabet_len)) / np.log(2)) * ppm
          - (np.log(background) * background / np.log(2))[None, :])
    return np.sum(ic, axis=1)


@pytest.mark.parametrize("key,name", sorted(RECORDED))
def test_ic_scale_recorded_values(key, name):
    pwm = ref_pwms(key)[name]
    ppm = pwm / np.sum(pwm, axis=-1, keepdims=True)
    total, ic0 = RECORDED[(key, name)]
    scaled = viz_sequence.ic_scale(ppm, background=[.25] * 4)
    assert scaled.shape == (20, 4)
    assert float(scaled.sum()) == pytest.approx(total, rel=1e-12)
    ic = viz_sequence.compute_per_position_ic(ppm, [.25] * 4, 0.001)
    assert float(ic[0]) == pytest.approx(ic0, rel=1e-12)
    # auto_shift_detect's wrapper (renormalise + uniform background) is the same computation
    np.testing.assert_array_equal(auto_shift_detect.ic_scale(pwm), scaled)


@pytest.mark.parametrize("key", ["atac_ref_motifs", "dnase_ref_motifs"])
def test_ic_scale_matches_modisco_formula_on_all_ref_motifs(key):
    for name, pwm in ref_pwms(key).items():
        for bg in ([.25] * 4, np.array([.3, .2, .2, .3])):
            want = pwm * modisco_0516_ic(pwm, np.asarray(bg))[:, None]
            np.testing.assert_array_equal(viz_sequence.ic_scale(pwm, background=bg), want, err_msg=name)


def test_background_list_and_array_agree():
    pwm = ref_pwms("atac_ref_motifs")["GSE101074_naive_hESC_ATAC_plus"]
    np.testing.assert_array_equal(viz_sequence.ic_scale(pwm, [.25] * 4),
                                  viz_sequence.ic_scale(pwm, np.full(4, .25)))


def test_row_renormalisation_semantics(capsys):
    # rows off by more than 1e-5 are renormalised for the IC, but the *input* pwm is what gets scaled
    x = np.array([[0.5, 0.5, 0.5, 0.5], [0.1, 0.2, 0.3, 0.4]])
    out = viz_sequence.ic_scale(x, [.3, .2, .2, .3])
    np.testing.assert_allclose(out, [[-0.014524702772665654] * 4,
                                     [0.012450780062474055, 0.02490156012494811,
                                      0.03735234018742216, 0.04980312024989622]], rtol=1e-12)
    assert "Will renormalize" in capsys.readouterr().out
    # within tolerance: no renormalisation, no warning
    y = np.array([[0.25, 0.25, 0.25, 0.25 + 5e-6]])
    ic = viz_sequence.compute_per_position_ic(y, [.25] * 4, 0.001)
    assert capsys.readouterr().out == ""
    raw = (np.log((y + 0.001) / 1.004) / np.log(2)) * y - (np.log(.25) * .25 / np.log(2))
    assert float(ic[0]) == float(raw.sum())


def test_ic_scale_checks_alphabet_axis():
    with pytest.raises(AssertionError):
        viz_sequence.ic_scale(np.full((4, 20), .25), [.25] * 3)


def test_plot_weights_given_ax_modisco_defaults():
    arr = np.zeros((10, 4))
    arr[2, 0] = 1.0
    arr[5, 3] = -0.5
    fig, ax = plt.subplots()
    viz_sequence.plot_weights_given_ax(ax=ax, array=arr)
    assert ax.get_xlim() == (-1.0, 11.0)
    np.testing.assert_allclose(ax.get_ylim(), (-0.5 - 0.2, 1.0 + 0.2))
    np.testing.assert_array_equal(ax.get_xticks(), np.arange(0.0, 11.0, 1.0))
    assert len(ax.patches) > 0
    plt.close(fig)

    # (4, L) input is transposed; (1, L, 4) is squeezed; highlight boxes and ylim are honoured
    fig, ax = plt.subplots()
    viz_sequence.plot_weights_given_ax(ax=ax, array=arr.T[None].transpose(0, 2, 1), subticks_frequency=5,
                                       highlight={"red": [(1, 3)]}, ylim=(-2.0, 2.0))
    np.testing.assert_allclose(ax.get_ylim(), (-2.4, 2.4))
    plt.close(fig)
    fig, ax = plt.subplots()
    viz_sequence.plot_weights_given_ax(ax=ax, array=arr.T)
    assert ax.get_xlim() == (-1.0, 11.0)
    plt.close(fig)


def test_plot_weights_returns_figure():
    fig = viz_sequence.plot_weights(np.eye(4)[np.arange(12) % 4] * 0.5)
    assert fig.get_size_inches().tolist() == [20.0, 2.0]
    plt.close(fig)
