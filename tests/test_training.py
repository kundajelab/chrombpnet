"""train.main end to end on a tiny synthetic genome / bigwig (bias model, then ChromBPNet with that bias model)."""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND=jax before keras is imported)
import jax
import keras
import matplotlib

matplotlib.use("Agg")

import chrombpnet.training.models.bpnet_model as bpnet_model
import chrombpnet.training.models.chrombpnet_with_bias_model as chrombpnet_with_bias_model
import chrombpnet.training.train as train
from chrombpnet.training.data_generators import initializers
from chrombpnet.training.data_generators.batchgen_generator import ChromBPNetBatchGenerator
from chrombpnet.training.utils import model_io
from chrombpnet.training.utils.callbacks import LossHistory
from chrombpnet.training.utils.losses import multinomial_nll

CHROMS = ("chr1", "chr2", "chr3")
CHROM_LEN = 6000
INPUTLEN, OUTPUTLEN, MAX_JITTER = 514, 200, 16
# parsers.py default
TRACKABLES = ["logcount_predictions_loss", "loss", "logits_profile_predictions_loss", "val_logcount_predictions_loss",
              "val_loss", "val_logits_profile_predictions_loss"]


def write_bed(path, chrom_summits):
    rows = [(c, s - 100, s + 100, ".", 0, ".", 0.0, 0.0, 0.0, 100) for c, s in chrom_summits]
    pd.DataFrame(rows).to_csv(path, sep="\t", header=False, index=False)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    import pyBigWig
    d = tmp_path_factory.mktemp("data")
    rng = np.random.RandomState(0)
    seqs = {c: "".join(rng.choice(list("ACGT"), CHROM_LEN)) for c in CHROMS}
    with open(d / "genome.fa", "w") as f:
        for c in CHROMS:
            f.write(">{}\n".format(c))
            f.writelines(seqs[c][i:i + 60] + "\n" for i in range(0, CHROM_LEN, 60))
    bw = pyBigWig.open(str(d / "signal.bw"), "w")
    bw.addHeader([(c, CHROM_LEN) for c in CHROMS], maxZooms=0)
    for c in CHROMS:
        gc = np.array([b in "GC" for b in seqs[c]], dtype=float)
        rate = 0.2 + 2 * np.convolve(gc, np.ones(9) / 9, mode="same")
        bw.addEntries(c, 0, values=rng.poisson(rate).astype(float).tolist(), span=1, step=1)
    bw.close()
    edge = INPUTLEN // 2 + MAX_JITTER + 10
    write_bed(d / "peaks.bed", [(c, s) for c in CHROMS for s in range(edge, CHROM_LEN - edge, 180)])
    write_bed(d / "nonpeaks.bed", [(c, s) for c in CHROMS for s in range(edge + 90, CHROM_LEN - edge, 180)])
    with open(d / "fold.json", "w") as f:
        json.dump({"train": ["chr1"], "valid": ["chr2"], "test": ["chr3"]}, f)
    return d


def write_params(path, data, **extra):
    params = {"counts_loss_weight": "10", "filters": "8", "n_dil_layers": "2", "inputlen": str(INPUTLEN),
              "outputlen": str(OUTPUTLEN), "negative_sampling_ratio": "0.5", "max_jitter": str(MAX_JITTER),
              "chr_fold_path": str(data / "fold.json")}
    params.update(extra)
    with open(path, "w") as f:
        f.write("\n".join("{}\t{}".format(k, v) for k, v in params.items()) + "\n")
    return str(path)


def make_args(data, out, architecture, peaks, nonpeaks, **extra):
    # the attributes the chrombpnet 1.x CLI passes to train.main; new options are optional
    args = argparse.Namespace(
        genome=str(data / "genome.fa"), bigwig=str(data / "signal.bw"), peaks=peaks, nonpeaks=nonpeaks,
        output_prefix=str(out), chr_fold_path=str(data / "fold.json"), epochs=2, early_stop=5, batch_size=8,
        learning_rate=1e-3, trackables=list(TRACKABLES), seed=1234, architecture_from_file=architecture,
        params=None, cmd="train")
    for k, v in extra.items():
        setattr(args, k, v)
    return args


def train_bias(data, out, **extra):
    args = make_args(data, out, bpnet_model.__file__, "None", str(data / "nonpeaks.bed"), **extra)
    args.params = write_params(str(out) + ".params.tsv", data)
    train.main(args)
    return args


def validation_loss(args, model_h5):
    """Loss of a saved model on the validation generator, as fit computes val_loss."""
    params = train.get_model_param_dict(args)
    model = model_io.load_model(model_h5)
    model.compile(loss=[multinomial_nll, "mse"], loss_weights=[1, float(params["counts_loss_weight"])])
    valid = initializers.initialize_generators(args, "valid", params, return_coords=False)
    return model.evaluate(valid, verbose=0, return_dict=True)["loss"]


@pytest.fixture(scope="module")
def bias_run(data, tmp_path_factory):
    out = tmp_path_factory.mktemp("bias") / "bias"
    return train_bias(data, out)


def test_bias_training_outputs(bias_run):
    prefix = bias_run.output_prefix
    for suffix in (".h5", ".log", ".log.batch", ".args.json"):
        assert os.path.exists(prefix + suffix), suffix
    log = pd.read_csv(prefix + ".log")
    assert list(log["epoch"]) == [0, 1]
    assert {"loss", "val_loss", "logits_profile_predictions_loss", "logcount_predictions_loss",
            "val_logits_profile_predictions_loss", "val_logcount_predictions_loss"} <= set(log.columns)
    assert np.isfinite(log.drop(columns="epoch").to_numpy()).all()
    batch = pd.read_csv(prefix + ".log.batch", sep="\t")
    assert list(batch.columns) == ["Epoch", "Batch"] + TRACKABLES
    assert np.isfinite(batch["loss"]).all() and set(batch["Epoch"]) == {0, 1}

    saved = json.load(open(prefix + ".args.json"))
    assert saved["seed"] == 1234 and saved["trackables"] == TRACKABLES
    assert saved["runtime"]["keras_backend"] == "jax" and saved["runtime"]["dtype_policy"] == "float32"
    # the architecture file is loaded under a real module name
    assert "chrombpnet_architecture_bpnet_model" in sys.modules


def test_saved_model_is_the_best_epoch(bias_run):
    # EarlyStopping(restore_best_weights) restores the best epoch at the end even when --epochs is reached
    log = pd.read_csv(bias_run.output_prefix + ".log")
    assert validation_loss(bias_run, bias_run.output_prefix + ".h5") == pytest.approx(log["val_loss"].min(), rel=1e-4)


def test_ema_training_saves_ema_weights(data, bias_run, tmp_path):
    args = train_bias(data, tmp_path / "bias_ema", ema=True, ema_momentum=0.9)
    log = pd.read_csv(args.output_prefix + ".log")
    # val_loss is computed on the EMA weights and the saved model holds them
    assert validation_loss(args, args.output_prefix + ".h5") == pytest.approx(log["val_loss"].min(), rel=1e-4)
    ema = model_io.load_model(args.output_prefix + ".h5").get_weights()
    plain = model_io.load_model(bias_run.output_prefix + ".h5").get_weights()
    assert not all(np.allclose(a, b) for a, b in zip(ema, plain))


def test_chrombpnet_training_with_bias(data, bias_run, tmp_path):
    out = tmp_path / "chrombpnet"
    args = make_args(data, out, chrombpnet_with_bias_model.__file__, str(data / "peaks.bed"),
                     str(data / "nonpeaks.bed"))
    args.params = write_params(str(out) + ".params.tsv", data, bias_model_path=bias_run.output_prefix + ".h5")
    train.main(args)
    prefix = args.output_prefix
    for suffix in (".h5", "_nobias.h5", ".log", ".log.batch", ".args.json"):
        assert os.path.exists(prefix + suffix), suffix
    full = model_io.load_model(prefix + ".h5")
    nobias = model_io.load_model(prefix + "_nobias.h5")
    bias = model_io.load_model(bias_run.output_prefix + ".h5")
    assert nobias.name == "model_wo_bias"
    x = initializers.initialize_generators(args, "valid", train.get_model_param_dict(args), return_coords=False)[0][0]
    logits, logcounts = full.predict(x, verbose=0)
    nb_logits, nb_logcounts = nobias.predict(x, verbose=0)
    b_logits, b_logcounts = bias.predict(x, verbose=0)
    # the bias model was frozen: full = nobias + bias
    np.testing.assert_allclose(logits, nb_logits + b_logits, rtol=1e-5, atol=1e-5)
    np.testing.assert_allclose(logcounts, np.logaddexp(nb_logcounts, b_logcounts), rtol=1e-5, atol=1e-5)
    log = pd.read_csv(prefix + ".log")
    assert validation_loss(args, prefix + ".h5") == pytest.approx(log["val_loss"].min(), rel=1e-4)


def test_opt_ins_muon_cosine_bf16(data, tmp_path):
    args = train_bias(data, tmp_path / "bias_opt", optimizer="muon", lr_schedule="cosine", precision="bf16",
                      device="auto")
    assert keras.config.dtype_policy().name == "float32"  # restored for the next pipeline steps
    saved = json.load(open(args.output_prefix + ".args.json"))
    assert saved["optimizer"] == "muon" and saved["precision"] == "bf16"
    # the cosine schedule spans --epochs full passes over the training batches
    train_batches = len(initializers.initialize_generators(args, "train", train.get_model_param_dict(args), False))
    assert saved["total_steps"] == args.epochs * train_batches > 0
    model = model_io.load_model(args.output_prefix + ".h5")
    # saved in float32 so that predict / interpret run it in full precision
    assert {layer.dtype_policy.name for layer in model.layers} == {"float32"}
    assert np.isfinite(pd.read_csv(args.output_prefix + ".log")["val_loss"]).all()


def test_precision_highest_is_restored(data, tmp_path):
    assert jax.config.jax_default_matmul_precision is None
    args = train_bias(data, tmp_path / "bias_highest", precision="highest", epochs=1)
    saved = json.load(open(args.output_prefix + ".args.json"))
    assert saved["runtime"]["jax_default_matmul_precision"] == "highest"
    # restored for the next pipeline steps (predict, interpret) and later runs in the same process
    assert jax.config.jax_default_matmul_precision is None


class Interrupted(Exception):
    pass


def test_interrupted_bf16_run_leaves_a_float32_checkpoint(data, tmp_path, monkeypatch):
    class InterruptAfterFirstEpoch(LossHistory):
        # runs after the ModelCheckpoint callback, as a time limit hit during epoch 2 would
        def on_epoch_end(self, epoch, logs=None):
            super().on_epoch_end(epoch, logs)
            assert self.model.get_layer("bpnet_1conv").compute_dtype == "bfloat16"  # still training in bf16
            raise Interrupted()

    monkeypatch.setattr(train.callbacks, "LossHistory", InterruptAfterFirstEpoch)
    out = tmp_path / "bias_bf16_interrupted"
    with pytest.raises(Interrupted):
        train_bias(data, out, precision="bf16")
    assert keras.config.dtype_policy().name == "float32"
    model = model_io.load_model(str(out) + ".h5")  # the epoch-1 checkpoint, not a final save
    policies = {layer.dtype_policy.name for layer in model._flatten_layers(include_self=True, recursive=True)}
    assert policies == {"float32"}
    assert model.get_layer("bpnet_1conv").compute_dtype == "float32"


def test_generator_contract(data):
    args = make_args(data, "unused", bpnet_model.__file__, str(data / "peaks.bed"), str(data / "nonpeaks.bed"),
                     inputlen=INPUTLEN, outputlen=OUTPUTLEN)
    params = {"inputlen": str(INPUTLEN), "outputlen": str(OUTPUTLEN), "negative_sampling_ratio": "0.5",
              "max_jitter": str(MAX_JITTER)}
    gen = initializers.initialize_generators(args, "train", params, return_coords=False)
    assert isinstance(gen, keras.utils.PyDataset) and isinstance(gen, ChromBPNetBatchGenerator)
    assert gen.workers == 1 and gen.use_multiprocessing is False
    x, y = gen[0]
    assert isinstance(y, tuple) and len(y) == 2
    assert x.dtype == np.int8 and x.shape == (8, INPUTLEN, 4)
    assert y[0].shape == (8, OUTPUTLEN) and y[1].shape == (8, 1)
    np.testing.assert_allclose(y[1], np.log(1 + y[0].sum(-1, keepdims=True)))
    assert len(gen) == int(np.ceil(gen.cur_seqs.shape[0] / 8))

    test = initializers.initialize_generators(args, "test", None, return_coords=True)
    x, y, coords = test[0]
    assert isinstance(y, tuple) and coords.shape == (8, 4)

    # optional thread workers (never processes) feeding fit
    threaded = initializers.initialize_generators(args, "train", params, return_coords=False, workers=2)
    assert threaded.workers == 2 and threaded.use_multiprocessing is False
    model = bpnet_model.getModelGivenModelOptionsAndWeightInits(
        argparse.Namespace(seed=1, learning_rate=1e-3),
        dict(params, filters="8", n_dil_layers="2", counts_loss_weight="10"))
    history = model.fit(threaded, epochs=2, verbose=0)
    assert np.isfinite(history.history["loss"]).all()
    with pytest.raises(ValueError):
        initializers.fetch_data_and_model_params_based_on_mode("predict", args, params, None, None)


def test_json_safe_args():
    args = argparse.Namespace(a=1, b="x", c=[1, "y"], d=None, e=True, f=object(), g=(1.5,), h={"k": 1})
    out = train.json_safe_args(args)
    json.dumps(out)
    assert out["a"] == 1 and out["c"] == [1, "y"] and out["g"] == [1.5] and out["d"] is None
    assert isinstance(out["f"], str) and isinstance(out["h"], str)


def test_loss_history_without_logs(tmp_path):
    history = LossHistory(str(tmp_path / "log.batch"), ["loss"])
    history.on_train_begin()
    history.on_epoch_begin(0)
    history.on_batch_end(0)
    history.on_batch_end(1, {"loss": 1.5})
    history.on_epoch_end(0)
    history.on_train_end()
    assert open(tmp_path / "log.batch").read().splitlines() == ["Epoch\tBatch\tloss", "0\t0\tNone", "0\t1\t1.5"]


def test_density_scatter_with_nans():
    from matplotlib import pyplot as plt
    from chrombpnet.training.utils.metrics_utils import density_scatter
    rng = np.random.RandomState(0)
    x, y = rng.normal(size=200), rng.normal(size=200)
    x[3] = np.nan
    y[[10, 20]] = np.nan  # different nan counts: a ragged index pair, which NumPy 2 rejects in np.isin
    plt.figure()
    density_scatter(x, y, xlab="x", ylab="y")
    plt.close("all")
