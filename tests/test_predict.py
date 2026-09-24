"""predict.main / predict_to_bigwig.main / reformat_chrombpnet_h5 end to end on tiny data (CPU)."""
import argparse
import json
import sys
import types

import h5py
import numpy as np
import pandas as pd
import pyBigWig
import pyfaidx
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND before keras is imported)
import keras
from scipy.special import logsumexp, softmax

from chrombpnet.training.utils import data_utils, model_io

INPUTLEN, OUTPUTLEN = 512, 256  # conv1 (21) + 2 dilated layers (2, 4) + profile conv (75): 512 -> 406 -> crop 256
CHROMS = [("chr1", 20000), ("chr2", 20000), ("chr3", 20000)]
PEAK_CENTERS = range(1500, 19000, 1500)
NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


def tiny_bpnet(filters=8, name=None, seed=0):
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
    counts = L.Dense(1, name="logcount_predictions")(L.GlobalAveragePooling1D(name="gap")(x))
    return keras.Model(inputs=inp, outputs=[prof, counts], name=name)


def write_dataset(d, rng=None):
    """Random genome, a Poisson bigWig (no zero bins) with extra reads at the summits, peaks/nonpeaks and a fold."""
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


def float_leaves(tree):
    if isinstance(tree, dict):
        return [x for v in tree.values() for x in float_leaves(v)]
    return [tree]


@pytest.fixture(scope="module")
def dataset(tmp_path_factory):
    d = tmp_path_factory.mktemp("predict")
    write_dataset(d)
    bias, nobias = tiny_bpnet(name="bias_model", seed=0), tiny_bpnet(filters=16, name="nobias", seed=1)
    bias.save(str(d / "bias.h5"))
    nobias.save(str(d / "chrombpnet_nobias.h5"))
    from chrombpnet.helpers.postprocessing import reformat_chrombpnet_h5
    reformat_chrombpnet_h5.main(str(d / "chrombpnet_nobias.h5"), str(d / "bias.h5"), str(d))
    (d / "chrombpnet_recompiled.h5").rename(d / "chrombpnet.h5")
    return d


@pytest.fixture
def predict_module(request, monkeypatch):
    """chrombpnet.training.predict, importing its data generator even before that generator is ported."""
    import chrombpnet.training.predict as predict
    try:
        import chrombpnet.training.data_generators.initializers  # noqa: F401
    except ImportError as e:
        if getattr(e, "name", None) != "tensorflow":
            raise
        # TEMPORARY (until batchgen_generator's Keras 3 port is merged): its `from tensorflow import keras` gets
        # Keras 3, whose keras.utils.Sequence is PyDataset; the generator is only indexed directly here
        stub = types.ModuleType("tensorflow")
        stub.keras = keras
        monkeypatch.setitem(sys.modules, "tensorflow", stub)
        generator_modules = ["chrombpnet.training.data_generators.batchgen_generator",
                             "chrombpnet.training.data_generators.initializers"]
        import chrombpnet.training.data_generators.initializers  # noqa: F401,F811
        request.addfinalizer(lambda: [sys.modules.pop(name, None) for name in generator_modules])
    return predict


def test_write_predictions_h5py_layout(tmp_path):
    from chrombpnet.training import predict
    profs = softmax(np.random.RandomState(0).randn(3, 8).astype(np.float32), axis=1)
    coords = np.array([["chr1", "100", "f", "1"], ["chr2", "250", "f", "0"], ["chrUn_x", "7", "r", "1"]])
    predict.write_predictions_h5py(str(tmp_path / "x"), profs, np.array([1.5, 2.5, 3.5], np.float32), coords)
    with h5py.File(tmp_path / "x_predictions.h5", "r") as f:
        assert sorted(f) == ["coords", "predictions"]
        assert sorted(f["coords"]) == ["coords_center", "coords_chrom", "coords_peak"]
        assert sorted(f["predictions"]) == ["logcounts", "profs"]
        assert h5py.check_string_dtype(f["coords/coords_chrom"].dtype).encoding == "utf-8"
        assert f["coords/coords_chrom"].asstr()[:].tolist() == ["chr1", "chr2", "chrUn_x"]
        assert f["coords/coords_center"].dtype == np.int64 and f["coords/coords_center"][:].tolist() == [100, 250, 7]
        assert f["coords/coords_peak"].dtype == np.int64 and f["coords/coords_peak"][:].tolist() == [1, 0, 1]
        assert f["predictions/profs"].dtype == np.float64 and f["predictions/profs"].shape == (3, 8)
        assert f["predictions/logcounts"].dtype == np.float64 and f["predictions/logcounts"].shape == (3,)
        assert all(f[k].compression == "gzip" for k in ["coords/coords_chrom", "predictions/profs"])


class BatchedArrays:
    def __init__(self, x, cts, batch_size):
        self.x, self.cts, self.batch_size = x, cts, batch_size
        self.coords = np.array([["chr1", str(i), "f", "1"] for i in range(len(x))])

    def __len__(self):
        return -(-len(self.x) // self.batch_size)

    def __getitem__(self, i):
        s = slice(i * self.batch_size, (i + 1) * self.batch_size)
        return self.x[s], [self.cts[s], np.log1p(self.cts[s].sum(-1, keepdims=True))], self.coords[s]


class ShapeSpy:
    def __init__(self, model):
        self.model, self.shapes = model, []

    def predict_on_batch(self, x):
        self.shapes.append(x.shape)
        return self.model.predict_on_batch(x)


@pytest.mark.parametrize("n,batch_size", [(10, 4), (3, 4), (8, 4)])
def test_predict_on_batch_wrapper_pads_short_batches(n, batch_size):
    from chrombpnet.training import predict
    rng = np.random.RandomState(0)
    x = np.eye(4, dtype=np.int8)[rng.randint(0, 4, (n, INPUTLEN))]
    cts = rng.poisson(2, (n, OUTPUTLEN)).astype(np.float64)
    model = tiny_bpnet()
    spy = ShapeSpy(model)
    true_counts, probs, true_sum, pred_sum, coords = predict.predict_on_batch_wrapper(spy, BatchedArrays(x, cts,
                                                                                                        batch_size))
    assert set(spy.shapes) == {(batch_size, INPUTLEN, 4)}  # one input shape -> one XLA compilation
    logits, logcounts = model.predict(x, verbose=0)
    assert probs.shape == (n, OUTPUTLEN) and pred_sum.shape == (n,) and coords.shape == (n, 4)
    np.testing.assert_allclose(probs, softmax(logits, axis=1), rtol=1e-5, atol=1e-7)
    np.testing.assert_allclose(pred_sum, logcounts[:, 0], rtol=1e-5, atol=1e-6)
    np.testing.assert_array_equal(true_counts, cts)
    np.testing.assert_allclose(true_sum, np.log1p(cts.sum(-1)))


def test_predict_main(dataset, tmp_path, predict_module):
    args = argparse.Namespace(model_h5=str(dataset / "bias.h5"), peaks=str(dataset / "peaks.bed"),
                              nonpeaks=str(dataset / "nonpeaks.bed"), output_prefix=str(tmp_path / "bias"),
                              batch_size=5, genome=str(dataset / "genome.fa"), bigwig=str(dataset / "signal.bw"),
                              chr_fold_path=str(dataset / "fold.json"), inputlen=INPUTLEN, outputlen=OUTPUTLEN)
    predict_module.main(args)

    peaks = pd.read_csv(dataset / "peaks.bed", sep="\t", names=NARROWPEAK_SCHEMA)
    nonpeaks = pd.read_csv(dataset / "nonpeaks.bed", sep="\t", names=NARROWPEAK_SCHEMA)
    test = pd.concat([peaks[peaks.chr == "chr3"], nonpeaks[nonpeaks.chr == "chr3"]])  # generator order
    with pyfaidx.Fasta(str(dataset / "genome.fa")) as genome:
        seqs = data_utils.get_seq(test, genome, INPUTLEN)
    logits, logcounts = model_io.load_model(dataset / "bias.h5").predict(seqs, verbose=0)
    with h5py.File(tmp_path / "bias_predictions.h5", "r") as f:
        assert f["coords/coords_chrom"].asstr()[:].tolist() == ["chr3"] * len(test)
        assert f["coords/coords_center"][:].tolist() == (test.start + test.summit).tolist()
        assert f["coords/coords_peak"][:].tolist() == [1] * 12 + [0] * 12
        profs = f["predictions/profs"][:]
        assert profs.dtype == np.float64 and profs.shape == (len(test), OUTPUTLEN)
        np.testing.assert_allclose(profs.sum(1), 1.0, rtol=1e-5)
        np.testing.assert_allclose(profs, softmax(logits, axis=1), rtol=1e-4, atol=1e-7)
        np.testing.assert_allclose(f["predictions/logcounts"][:], logcounts[:, 0], rtol=1e-5, atol=1e-6)

    metrics = json.load(open(tmp_path / "bias_metrics.json"))
    assert sorted(metrics["counts_metrics"]) == ["nonpeaks", "peaks", "peaks_and_nonpeaks"]
    assert sorted(metrics["counts_metrics"]["peaks"]) == ["mse", "pearsonr", "spearmanr"]
    assert sorted(metrics["profile_metrics"]["peaks"]) == ["median_jsd", "median_norm_jsd"]
    leaves = float_leaves(metrics)
    assert len(leaves) == 15 and all(type(v) is float and np.isfinite(v) for v in leaves)
    with pyBigWig.open(str(dataset / "signal.bw")) as bw:
        labels = np.log1p(data_utils.get_cts(test, bw, OUTPUTLEN).sum(-1))
    np.testing.assert_allclose(metrics["counts_metrics"]["peaks_and_nonpeaks"]["mse"],
                               np.mean((labels - logcounts[:, 0]) ** 2), rtol=1e-4)
    for name in ["_peaks_and_nonpeaks", "_only_peaks", "_only_nonpeaks"]:
        assert (tmp_path / ("bias" + name + ".counts_pearsonr.png")).exists()
        assert (tmp_path / ("bias" + name + ".profile_jsd.png")).exists()


def test_predict_main_peaks_only(dataset, tmp_path, predict_module):
    args = argparse.Namespace(model_h5=str(dataset / "chrombpnet.h5"), peaks=str(dataset / "peaks.bed"),
                              nonpeaks="None", output_prefix=str(tmp_path / "chrombpnet"), batch_size=64,
                              genome=str(dataset / "genome.fa"), bigwig=str(dataset / "signal.bw"),
                              chr_fold_path=str(dataset / "fold.json"), inputlen=INPUTLEN, outputlen=OUTPUTLEN)
    predict_module.main(args)
    metrics = json.load(open(tmp_path / "chrombpnet_metrics.json"))
    assert list(metrics["counts_metrics"]) == ["peaks"]
    assert all(type(v) is float for v in float_leaves(metrics))


def test_reformat_chrombpnet_h5(dataset):
    full = model_io.load_model(dataset / "chrombpnet.h5")
    bias = model_io.load_model(dataset / "bias.h5")
    nobias = model_io.load_model(dataset / "chrombpnet_nobias.h5")
    assert sorted(layer.name for layer in full.layers) == ["bias_model", "concatenate", "logcount_predictions",
                                                           "logits_profile_predictions", "model_wo_bias", "sequence"]
    assert type(full.get_layer("logcount_predictions")).__name__ == "LogSumExp"
    seqs = np.eye(4, dtype=np.int8)[np.random.RandomState(3).randint(0, 4, (6, INPUTLEN))]
    (fl, fc), (bl, bc), (nl, nc) = [m.predict(seqs, verbose=0) for m in (full, bias, nobias)]
    np.testing.assert_allclose(fl, nl + bl, rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(fc, logsumexp(np.concatenate([nc, bc], axis=1), axis=1, keepdims=True), rtol=1e-5,
                               atol=1e-5)
    assert not (dataset / "chrombpnet_recompiled").exists()  # no SavedModel export any more


def test_predict_to_bigwig(dataset, tmp_path):
    from chrombpnet.evaluation.make_bigwigs import predict_to_bigwig
    prefix = str(tmp_path / "pred")
    args = argparse.Namespace(bias_model=str(dataset / "bias.h5"), chrombpnet_model=str(dataset / "chrombpnet.h5"),
                              chrombpnet_model_nb=str(dataset / "chrombpnet_nobias.h5"),
                              regions=str(dataset / "peaks.bed"), genome=str(dataset / "genome.fa"),
                              chrom_sizes=str(dataset / "chrom.sizes"), output_prefix=prefix,
                              output_prefix_stats=str(tmp_path / "stats.txt"), batch_size=7, tqdm=0, debug_chr=None,
                              bigwig=str(dataset / "signal.bw"))
    predict_to_bigwig.main(args)

    regions = pd.read_csv(dataset / "peaks.bed", sep="\t", names=NARROWPEAK_SCHEMA)
    with pyfaidx.Fasta(str(dataset / "genome.fa")) as genome:
        seqs = data_utils.get_seq(regions, genome, INPUTLEN)
    for name, path in [("bias", "bias.h5"), ("chrombpnet", "chrombpnet.h5"),
                       ("chrombpnet_nobias", "chrombpnet_nobias.h5")]:
        logits, logcounts = model_io.load_model(dataset / path).predict(seqs, verbose=0)
        expected = softmax(logits, axis=1) * np.exp(logcounts)
        with pyBigWig.open(prefix + "_{}.bw".format(name)) as bw:
            assert bw.chroms() == dict(CHROMS)
            for i in (0, len(regions) - 1):
                r = regions.iloc[i]
                mid = r.start + r.summit
                got = np.array(bw.values(r.chr, mid - OUTPUTLEN // 2, mid + OUTPUTLEN // 2))
                np.testing.assert_allclose(got, expected[i], rtol=1e-4)
        assert len(pd.read_csv(prefix + "_{}_preds.bed".format(name), sep="\t", header=None)) == len(regions)
        with h5py.File(prefix + "_{}_predictions.h5".format(name), "r") as f:
            assert sorted(f["coords"]) == ["coords_center", "coords_chrom"]
            assert f["predictions/profs"].dtype == np.float64
            np.testing.assert_allclose(f["predictions/logcounts"][:], logcounts[:, 0], rtol=1e-5, atol=1e-6)
        metrics = json.load(open(prefix + "_{}_metrics.json".format(name)))
        assert all(type(v) is float for v in float_leaves(metrics))
    assert open(tmp_path / "stats.txt").read().startswith("Min\t")


# ---------------------------------------------------------------- -d/--debug-chr on a 2-chromosome genome

TWO_CHROMS = [("chr1", 8000), ("chr2", 6000)]
# interleaved chromosomes; chr2:5900 runs past the end of chr2 and is dropped (regions_used)
TWO_CHROM_CENTERS = [("chr2", 1000), ("chr1", 1000), ("chr1", 3000), ("chr2", 5900), ("chr2", 3000), ("chr1", 5000),
                     ("chr2", 2000)]


@pytest.fixture(scope="module")
def two_chrom(tmp_path_factory):
    d = tmp_path_factory.mktemp("two_chrom")
    rng = np.random.RandomState(1)
    with open(d / "genome.fa", "w") as f:
        for c, n in TWO_CHROMS:
            seq = "".join(rng.choice(list("ACGT"), n))
            f.write(">{}\n".format(c))
            f.writelines(seq[i:i + 60] + "\n" for i in range(0, n, 60))
    (d / "chrom.sizes").write_text("".join("{}\t{}\n".format(c, n) for c, n in TWO_CHROMS))
    rows = [[c, m - 250, m + 250, ".", ".", ".", ".", ".", ".", 250] for c, m in TWO_CHROM_CENTERS]
    pd.DataFrame(rows).to_csv(d / "peaks.bed", sep="\t", header=False, index=False)
    return d


def _kept_regions(d, chroms):
    regions = pd.read_csv(d / "peaks.bed", sep="\t", names=NARROWPEAK_SCHEMA)
    mid = regions.start + regions.summit
    fits = (mid + INPUTLEN // 2) <= regions.chr.map(dict(TWO_CHROMS))
    return regions[regions.chr.isin(chroms) & fits].reset_index(drop=True)


@pytest.mark.parametrize("debug_chr", [["chr2"], ["chr1", "chr2"]])
def test_pred_bw_debug_chr(dataset, two_chrom, tmp_path, monkeypatch, debug_chr):
    import chrombpnet.CHROMBPNET as cli
    prefix = str(tmp_path / "pred")
    monkeypatch.setattr(sys, "argv", ["chrombpnet", "pred_bw", "-bm", str(dataset / "bias.h5"),
                                      "-cm", str(dataset / "chrombpnet.h5"),
                                      "-cmb", str(dataset / "chrombpnet_nobias.h5"), "-r", str(two_chrom / "peaks.bed"),
                                      "-g", str(two_chrom / "genome.fa"), "-c", str(two_chrom / "chrom.sizes"),
                                      "-op", prefix, "-t", "0", "-bs", "2", "-d"] + debug_chr)
    cli.main()

    kept = _kept_regions(two_chrom, debug_chr)
    assert len(kept) == (3 if debug_chr == ["chr2"] else 6)
    with pyfaidx.Fasta(str(two_chrom / "genome.fa")) as genome:
        seqs = data_utils.get_seq(kept, genome, INPUTLEN)
    for name, path in [("bias", "bias.h5"), ("chrombpnet", "chrombpnet.h5"),
                       ("chrombpnet_nobias", "chrombpnet_nobias.h5")]:
        preds_bed = pd.read_csv(prefix + "_{}_preds.bed".format(name), sep="\t", names=NARROWPEAK_SCHEMA)
        pd.testing.assert_frame_equal(preds_bed, kept)
        logits, logcounts = model_io.load_model(dataset / path).predict(seqs, verbose=0)
        expected = softmax(logits, axis=1) * np.exp(logcounts)
        with pyBigWig.open(prefix + "_{}.bw".format(name)) as bw:
            assert bw.chroms() == dict(TWO_CHROMS)
            for c, _ in TWO_CHROMS:
                if c not in debug_chr:
                    assert not bw.intervals(c)
            assert sum(len(bw.intervals(c) or ()) for c in debug_chr) == len(kept) * OUTPUTLEN
            for k, r in kept.iterrows():
                mid = r.start + r.summit
                got = np.array(bw.values(r.chr, mid - OUTPUTLEN // 2, mid + OUTPUTLEN // 2))
                np.testing.assert_allclose(got, expected[k], rtol=1e-4)


def test_contribs_bw_debug_chr(dataset, two_chrom, tmp_path, monkeypatch):
    import hdf5plugin  # noqa: F401  (Blosc filter)
    import chrombpnet.CHROMBPNET as cli
    prefix = str(tmp_path / "contribs")
    monkeypatch.setattr(sys, "argv", ["chrombpnet", "contribs_bw", "-m", str(dataset / "bias.h5"),
                                      "-r", str(two_chrom / "peaks.bed"), "-g", str(two_chrom / "genome.fa"),
                                      "-c", str(two_chrom / "chrom.sizes"), "-op", prefix, "-pc", "counts", "-t", "0",
                                      "--shap-batch-seqs", "2", "-d", "chr2"])
    cli.main()
    kept = _kept_regions(two_chrom, ["chr2"])
    interpreted = pd.read_csv(prefix + ".interpreted_regions.bed", sep="\t", names=NARROWPEAK_SCHEMA)
    pd.testing.assert_frame_equal(interpreted, kept)
    with h5py.File(prefix + ".counts_scores.h5", "r") as f:
        scores = f["projected_shap"]["seq"][:].sum(1)
    assert scores.shape == (len(kept), INPUTLEN)
    with pyBigWig.open(prefix + ".counts_scores.bw") as bw:
        assert not bw.intervals("chr1")
        for k, r in kept.iterrows():
            mid = r.start + r.summit
            np.testing.assert_array_equal(np.array(bw.values("chr2", mid - INPUTLEN // 2, mid + INPUTLEN // 2),
                                                   np.float32), scores[k].astype(np.float32))


@pytest.mark.parametrize("debug_chr", [None, "chr2", ["chr2"], ("chr1", "chr2")])
def test_write_bigwig_debug_chr(tmp_path, debug_chr):
    """write_bigwig subsets to one chromosome name (importance_hdf5_to_bigwig -d) or a list of them (contribs_bw)."""
    from chrombpnet.evaluation.make_bigwigs import bigwig_helper
    gs = [("chr1", 1000), ("chr2", 1000)]
    regions = [["chr2", 100, 110, 105], ["chr1", 200, 210, 205], ["chr2", 50, 60, 55]]
    data = np.arange(30, dtype=float).reshape(3, 10)
    bigwig_helper.write_bigwig(data, regions, gs, str(tmp_path / "x.bw"), debug_chr=debug_chr,
                               outstats_file=str(tmp_path / "stats.txt"))
    keep = {"chr1", "chr2"} if debug_chr is None else {debug_chr} if isinstance(debug_chr, str) else set(debug_chr)
    with pyBigWig.open(str(tmp_path / "x.bw")) as bw:
        for k, (c, s, e, _) in enumerate(regions):
            if c in keep:
                assert bw.values(c, s, e) == data[k].tolist()
            else:
                assert not bw.intervals(c)
    kept = np.concatenate([data[k] for k, r in enumerate(regions) if r[0] in keep])
    stats = dict(line.split("\t") for line in open(tmp_path / "stats.txt").read().splitlines())
    assert (float(stats["Min"]), float(stats["Max"])) == (kept.min(), kept.max())
