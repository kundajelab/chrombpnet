"""find_bias_hyperparams / find_chrombpnet_hyperparams on a synthetic genome + bigWig with a tiny Keras 3 bias model."""
import argparse
import json
import sys

import numpy as np
import pandas as pd
import pyBigWig
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND before keras is imported)
import keras

from chrombpnet.helpers.hyperparameters import find_bias_hyperparams, find_chrombpnet_hyperparams, param_utils
from chrombpnet.training.utils import model_io

INPUTLEN, OUTPUTLEN = 512, 256  # conv1 (21) + 2 dilated layers (2, 4) + profile conv (75): 512 -> 406 -> crop 256
CHROMS = [("chr1", 20000), ("chr2", 20000), ("chr3", 20000)]
PEAK_CENTERS = range(1500, 19000, 1500)


def tiny_bpnet(filters=8, count_head="logcount_predictions", name=None, seed=0):
    keras.utils.set_random_seed(seed)
    L = keras.layers
    inp = L.Input(shape=(INPUTLEN, 4), name="sequence")
    x = L.Conv1D(filters, 21, activation="relu", name="bpnet_1st_conv")(inp)
    for i in (1, 2):
        conv = L.Conv1D(filters, 3, activation="relu", dilation_rate=2 ** i, name="bpnet_{}conv".format(i))(x)
        x = L.Cropping1D((x.shape[1] - conv.shape[1]) // 2, name="bpnet_{}crop".format(i))(x)
        x = L.add([conv, x])
    prof = L.Conv1D(1, 75, name="prof_out_precrop")(x)
    prof = L.Cropping1D(prof.shape[1] // 2 - OUTPUTLEN // 2, name="logits_profile_predictions_preflatten")(prof)
    prof = L.Flatten(name="logits_profile_predictions")(prof)
    counts = L.Dense(1, name=count_head)(L.GlobalAveragePooling1D(name="gap")(x))
    return keras.Model(inputs=inp, outputs=[prof, counts], name=name)


def write_dataset(d, rng=None):
    """Random genome, a Poisson bigWig with extra reads around the peak summits, peaks/nonpeaks and a fold."""
    rng = rng or np.random.RandomState(0)
    with open(d / "genome.fa", "w") as f:
        for c, n in CHROMS:
            seq = "".join(rng.choice(list("ACGT"), n))
            f.write(">{}\n".format(c))
            f.writelines(seq[i:i + 60] + "\n" for i in range(0, n, 60))
    (d / "chrom.sizes").write_text("".join("{}\t{}\n".format(c, n) for c, n in CHROMS))
    bw = pyBigWig.open(str(d / "signal.bw"), "w")
    bw.addHeader(CHROMS)
    for c, n in CHROMS:
        vals = rng.poisson(1.0, n).astype(np.float64) + 1.0
        for m in PEAK_CENTERS:
            vals[m - 200:m + 200] += rng.poisson(rng.uniform(5, 15), 400)
        bw.addEntries(c, 0, values=vals, span=1, step=1)
    bw.close()
    peaks = [[c, m - 250, m + 250, ".", ".", ".", ".", ".", ".", 250] for c, _ in CHROMS for m in PEAK_CENTERS]
    nonpeaks = [[c, m + 500, m + 1000, ".", ".", ".", ".", ".", ".", 250] for c, _ in CHROMS for m in PEAK_CENTERS]
    pd.DataFrame(peaks).to_csv(d / "peaks.bed", sep="\t", header=False, index=False)
    pd.DataFrame(nonpeaks).to_csv(d / "nonpeaks.bed", sep="\t", header=False, index=False)
    json.dump({"train": ["chr1"], "valid": ["chr2"], "test": ["chr3"]}, open(d / "fold.json", "w"))


def read_tsv(path):
    return dict(line.rstrip("\n").split("\t") for line in open(path) if line.strip())


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    d = tmp_path_factory.mktemp("hyperparams")
    write_dataset(d)
    tiny_bpnet(name="bias_model").save(str(d / "bias.h5"))
    return d


def chrombpnet_args(d, out, **extra):
    base = dict(genome=str(d / "genome.fa"), bigwig=str(d / "signal.bw"), peaks=str(d / "peaks.bed"),
                nonpeaks=str(d / "nonpeaks.bed"), negative_sampling_ratio=0.5, outlier_threshold=0.9999, max_jitter=16,
                chr_fold_path=str(d / "fold.json"), inputlen=INPUTLEN, outputlen=OUTPUTLEN, filters=16,
                n_dilation_layers=2, bias_model_path=str(d / "bias.h5"), output_prefix=str(out) + "/", seed=1234)
    base.update(extra)
    return argparse.Namespace(**base)


def test_param_utils_does_not_import_keras():
    code = ("import sys; sys.modules['tensorflow'] = None\n"
            "import chrombpnet.helpers.hyperparameters.find_bias_hyperparams\n"
            "assert 'keras' not in sys.modules\n")
    import subprocess
    assert subprocess.run([sys.executable, "-c", code]).returncode == 0


def test_load_model_wrapper_alias(dataset):
    model = param_utils.load_model_wrapper(model_h5=str(dataset / "bias.h5"))
    assert model.input_shape == (None, INPUTLEN, 4)
    assert model.output_shape == [(None, OUTPUTLEN), (None, 1)]
    assert getattr(model, "optimizer", None) is None  # compile=False


def test_find_bias_hyperparams(dataset, tmp_path):
    args = argparse.Namespace(genome=str(dataset / "genome.fa"), bigwig=str(dataset / "signal.bw"),
                              peaks=str(dataset / "peaks.bed"), nonpeaks=str(dataset / "nonpeaks.bed"),
                              bias_threshold_factor=0.5, outlier_threshold=0.9999, max_jitter=0,
                              chr_fold_path=str(dataset / "fold.json"), inputlen=INPUTLEN, outputlen=OUTPUTLEN,
                              filters=8, n_dilation_layers=2, output_prefix=str(tmp_path) + "/")
    find_bias_hyperparams.main(args)
    params = read_tsv(tmp_path / "bias_model_params.tsv")
    assert list(params) == ["counts_loss_weight", "filters", "n_dil_layers", "inputlen", "outputlen", "max_jitter",
                            "chr_fold_path", "negative_sampling_ratio"]
    assert (params["filters"], params["inputlen"], params["outputlen"]) == ("8", str(INPUTLEN), str(OUTPUTLEN))
    assert float(params["counts_loss_weight"]) >= 1.0
    data = read_tsv(tmp_path / "bias_data_params.tsv")
    assert float(data["counts_sum_min_thresh"]) < float(data["counts_sum_max_thresh"])
    nonpeaks = pd.read_csv(tmp_path / "filtered.bias_nonpeaks.bed", sep="\t", header=None)
    # train/valid nonpeaks are far below the peak counts, so only the outliers (min/max, with ties) are removed;
    # test nonpeaks are never filtered
    assert (nonpeaks[0] == "chr3").sum() == len(PEAK_CENTERS)
    assert 2 * len(PEAK_CENTERS) - 4 <= (nonpeaks[0] != "chr3").sum() <= 2 * len(PEAK_CENTERS) - 2


def test_adjust_bias_model_logcounts():
    rng = np.random.RandomState(1)
    seqs = np.eye(4, dtype=np.int8)[rng.randint(0, 4, (24, INPUTLEN))]
    cts = rng.poisson(300, 24).astype(np.float64)
    for head in ("logcount_predictions", "logcounts"):
        model = tiny_bpnet(count_head=head)
        kernel, bias = [w.copy() for w in model.get_layer(head).get_weights()]
        before = model.predict(seqs, verbose=0)[1].ravel()
        scaled = find_chrombpnet_hyperparams.adjust_bias_model_logcounts(model, seqs, cts)
        after = scaled.predict(seqs, verbose=0)[1].ravel()
        new_kernel, new_bias = scaled.get_layer(head).get_weights()
        np.testing.assert_array_equal(new_kernel, kernel)
        assert new_bias.dtype == np.float32
        delta = np.mean(np.log1p(cts) - before)
        np.testing.assert_allclose(new_bias - bias, delta, rtol=1e-5)
        np.testing.assert_allclose(np.mean(after), np.mean(np.log1p(cts)), atol=1e-4)


def test_adjust_bias_model_logcounts_refuses_other_heads():
    model = tiny_bpnet(count_head="some_dense")
    with pytest.raises(AssertionError):
        find_chrombpnet_hyperparams.adjust_bias_model_logcounts(model, np.zeros((2, INPUTLEN, 4), np.int8),
                                                                np.ones(2))


def test_find_chrombpnet_hyperparams(dataset, tmp_path):
    runs = []
    for i in range(2):
        out = tmp_path / "run{}".format(i)
        out.mkdir()
        find_chrombpnet_hyperparams.main(chrombpnet_args(dataset, out))
        runs.append(out)
    out = runs[0]
    params = read_tsv(out / "chrombpnet_model_params.tsv")
    assert params["bias_model_path"] == str(out) + "/bias_model_scaled.h5"
    assert (params["filters"], params["n_dil_layers"], params["max_jitter"]) == ("16", "2", "16")
    assert (params["inputlen"], params["outputlen"], params["negative_sampling_ratio"]) == (str(INPUTLEN),
                                                                                            str(OUTPUTLEN), "0.5")
    assert float(params["counts_loss_weight"]) > 1.0
    # the negative subsample behind the outlier thresholds is seeded with --seed
    for name in ("chrombpnet_data_params.tsv", "filtered.peaks.bed", "filtered.nonpeaks.bed"):
        assert (runs[0] / name).read_text() == (runs[1] / name).read_text()

    original = model_io.load_model(dataset / "bias.h5")
    scaled = model_io.load_model(out / "bias_model_scaled.h5")
    assert scaled.input_shape == original.input_shape and scaled.output_shape == original.output_shape
    changed = [(layer.name, i) for layer in original.layers for i, (w0, w1) in
               enumerate(zip(layer.get_weights(), scaled.get_layer(layer.name).get_weights())) if not np.array_equal(w0, w1)]
    assert changed == [("logcount_predictions", 1)]  # only the count head's bias moves
    delta = (scaled.get_layer("logcount_predictions").get_weights()[1]
             - original.get_layer("logcount_predictions").get_weights()[1])
    seqs = np.eye(4, dtype=np.int8)[np.random.RandomState(2).randint(0, 4, (4, INPUTLEN))]
    shift = scaled.predict(seqs, verbose=0)[1] - original.predict(seqs, verbose=0)[1]
    np.testing.assert_allclose(shift, np.full((4, 1), delta[0]), atol=1e-5)
    assert abs(delta[0]) > 0.1


def test_find_chrombpnet_hyperparams_without_seed(dataset, tmp_path):
    # the standalone parser has no --seed: the subsample stays unseeded, as in chrombpnet 1.x
    args = chrombpnet_args(dataset, tmp_path)
    del args.seed
    find_chrombpnet_hyperparams.main(args)
    assert (tmp_path / "bias_model_scaled.h5").exists()
