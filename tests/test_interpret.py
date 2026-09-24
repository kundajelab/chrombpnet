"""interpret.main end to end on a synthetic genome with a tiny model: output files, h5 contract, references,
settings resolution, and the contribution-score bigwig."""
import argparse
import json
import os
import subprocess
import sys

import chrombpnet  # noqa: F401  (KERAS_BACKEND=jax)
import h5py
import hdf5plugin  # noqa: F401  (Blosc filter)
import keras
import numpy as np
import pandas as pd
import pyBigWig
import pytest
from keras import layers

import chrombpnet.evaluation.interpret.interpret as interpret
import chrombpnet.evaluation.make_bigwigs.importance_hdf5_to_bigwig as importance_hdf5_to_bigwig
from chrombpnet.evaluation.interpret import scores_io
from chrombpnet.evaluation.interpret.explainer import DeepLiftShap
from chrombpnet.training.utils.one_hot import dna_to_one_hot

L, OUTLEN = 256, 100
CHROMS = {"chr1": 4000, "chr2": 2000}
CENTERS = [("chr1", 600), ("chr1", 1200), ("chr2", 500), ("chr1", 1900), ("chr1", 3990), ("chr1", 2600),
           ("chr2", 1100), ("chr1", 3300)]
EDGE = 4  # chr1:3990 runs past the chromosome end and is dropped


def tiny_model(seed=0):
    keras.utils.set_random_seed(seed)
    inp = keras.Input((L, 4), name="sequence")
    x = layers.Conv1D(8, 21, activation="relu", name="bpnet_1st_conv")(inp)
    for i in (1, 2):
        conv_x = layers.Conv1D(8, 3, activation="relu", dilation_rate=2 ** i, name="bpnet_{}conv".format(i))(x)
        x = layers.Cropping1D((x.shape[1] - conv_x.shape[1]) // 2, name="bpnet_{}crop".format(i))(x)
        x = layers.add([conv_x, x])
    prof = layers.Conv1D(1, 75, name="prof_out_precrop")(x)
    prof = layers.Cropping1D(prof.shape[1] // 2 - OUTLEN // 2, name="logits_profile_predictions_preflatten")(prof)
    profile_out = layers.Flatten(name="logits_profile_predictions")(prof)
    count_out = layers.Dense(1, name="logcount_predictions")(layers.GlobalAveragePooling1D(name="gap")(x))
    model = keras.Model(inp, [profile_out, count_out])
    rng = np.random.RandomState(seed)
    model.set_weights([w if w.ndim > 1 else rng.normal(0, 0.05, w.shape).astype(np.float32)
                       for w in model.get_weights()])
    return model


@pytest.fixture(scope="module")
def workdir(tmp_path_factory):
    d = tmp_path_factory.mktemp("interpret")
    rng = np.random.RandomState(0)
    genome = {c: "".join(rng.choice(list("ACGT"), n)) for c, n in CHROMS.items()}
    genome["chr1"] = genome["chr1"][:1150] + "N" * 20 + genome["chr1"][1170:]  # N run inside a region
    with open(d / "genome.fa", "w") as f:
        for c, s in genome.items():
            f.write(">{}\n".format(c))
            for i in range(0, len(s), 60):
                f.write(s[i:i + 60] + "\n")
    with open(d / "chrom.sizes", "w") as f:
        for c, n in CHROMS.items():
            f.write("{}\t{}\n".format(c, n))
    rows = [[c, m - 50, m + 50, ".", 0, ".", 0, 0, 0, 50] for c, m in CENTERS]
    pd.DataFrame(rows).to_csv(d / "peaks.bed", sep="\t", header=False, index=False)
    tiny_model().save(str(d / "model.h5"))
    return {"dir": d, "genome": genome}


def ns(workdir, prefix, **kw):
    d = workdir["dir"]
    args = dict(genome=str(d / "genome.fa"), regions=str(d / "peaks.bed"), model_h5=str(d / "model.h5"),
                output_prefix=str(d / prefix), debug_chr=None, profile_or_counts=["counts", "profile"],
                batch_seqs=2)
    args.update(kw)
    return argparse.Namespace(**args)


def expected_seqs(workdir):
    g = workdir["genome"]
    return dna_to_one_hot([g[c][m - L // 2:m + L // 2] for i, (c, m) in enumerate(CENTERS) if i != EDGE])


def read(path):
    with h5py.File(path, "r") as f:
        return {k: f[k]["seq"][:] for k in scores_io.KEYS}


@pytest.fixture(scope="module")
def run_default(workdir):
    args = ns(workdir, "run", save_references=True)
    interpret.main(args)
    return args


def test_output_files_and_h5_contract(workdir, run_default):
    prefix = run_default.output_prefix
    for suffix in (".interpret.args.json", ".interpreted_regions.bed", ".counts_scores.h5", ".profile_scores.h5",
                   ".references.npz"):
        assert os.path.exists(prefix + suffix), suffix
    assert not any(n.endswith(".partial") for n in os.listdir(workdir["dir"]))
    bed = pd.read_csv(prefix + ".interpreted_regions.bed", sep="\t", header=None)
    assert len(bed) == len(CENTERS) - 1
    assert list(zip(bed[0], bed[1] + bed[9])) == [cm for i, cm in enumerate(CENTERS) if i != EDGE]
    x = expected_seqs(workdir)
    n = len(x)
    for head in ("counts", "profile"):
        path = "{}.{}_scores.h5".format(prefix, head)
        with h5py.File(path, "r") as f:
            assert set(f.keys()) == set(scores_io.KEYS)
            for key, dtype in scores_io.DTYPES.items():
                ds = f[key]["seq"]
                assert ds.shape == (n, 4, L) and ds.dtype == dtype
                assert "32001" in ds._filters  # Blosc, as deepdish/PyTables wrote it
        # what `modisco motifs -i` reads
        d = read(path)
        np.testing.assert_array_equal(d["raw"], np.transpose(x, (0, 2, 1)))
        np.testing.assert_array_equal(d["projected_shap"], d["raw"] * d["shap"])
        n_pos = x.sum(-1) == 0
        assert n_pos.sum() == 20 and np.all(np.transpose(d["projected_shap"], (0, 2, 1))[n_pos] == 0)
        assert np.abs(d["shap"]).max() > 0


def test_scores_match_explainer(workdir, run_default):
    x = expected_seqs(workdir)
    from chrombpnet.training.utils.model_io import load_model
    hyp = DeepLiftShap(load_model(run_default.model_h5), batch_seqs=2).explain(x, seed=1234)
    for head in ("counts", "profile"):
        d = read("{}.{}_scores.h5".format(run_default.output_prefix, head))
        np.testing.assert_array_equal(d["shap"], np.transpose(hyp[head], (0, 2, 1)).astype(np.float16))


def test_args_json_records_settings(run_default):
    with open(run_default.output_prefix + ".interpret.args.json") as f:
        rec = json.load(f)
    assert rec["precision"] == "auto" and rec["device"] == "auto" and rec["seed"] == 1234
    assert rec["precision_resolved"] == ("highest" if rec["jax_backend"] == "cpu" else "default")
    assert rec["batch_seqs"] == 2 and rec["num_shuffles"] == 20 and rec["profile_weighting"] == "chrombpnet"
    assert rec["jax_backend"] in ("cpu", "gpu") and rec["jax_devices"]
    assert rec["model_h5"] == run_default.model_h5


def test_saved_references_reproduce_run(workdir, run_default):
    prefix = run_default.output_prefix
    refs = np.load(prefix + ".references.npz")
    x = expected_seqs(workdir)
    np.testing.assert_array_equal(refs["x"], x)
    assert refs["refs"].shape == (len(x), 20, L, 4) and refs["refs"].dtype == np.int8
    args = ns(workdir, "rerun", references=prefix + ".references.npz", seed=99, batch_seqs=3)
    interpret.main(args)
    for head in ("counts", "profile"):
        a = read("{}.{}_scores.h5".format(prefix, head))
        b = read("{}.{}_scores.h5".format(args.output_prefix, head))
        np.testing.assert_array_equal(a["raw"], b["raw"])
        diff = np.abs(a["shap"].astype(np.float32) - b["shap"].astype(np.float32))
        assert diff.max() <= 2e-3 * np.abs(a["shap"].astype(np.float32)).max()


def test_region_order_and_seed(workdir, run_default):
    d = workdir["dir"]
    peaks = pd.read_csv(d / "peaks.bed", sep="\t", header=None)
    peaks.iloc[::-1].to_csv(d / "peaks_rev.bed", sep="\t", header=False, index=False)
    args = ns(workdir, "rev", regions=str(d / "peaks_rev.bed"), profile_or_counts=["counts"])
    interpret.main(args)
    assert not os.path.exists(args.output_prefix + ".profile_scores.h5")
    a = read(run_default.output_prefix + ".counts_scores.h5")
    b = read(args.output_prefix + ".counts_scores.h5")
    np.testing.assert_array_equal(b["raw"][::-1], a["raw"])
    np.testing.assert_allclose(b["shap"][::-1].astype(np.float32), a["shap"].astype(np.float32),
                               rtol=2e-3, atol=1e-3 * np.abs(a["shap"].astype(np.float32)).max())
    other = ns(workdir, "seed7", profile_or_counts=["counts"], shap_seed=7, seed=1234)
    interpret.main(other)
    c = read(other.output_prefix + ".counts_scores.h5")
    assert not np.array_equal(c["shap"], a["shap"])
    with open(other.output_prefix + ".interpret.args.json") as f:
        assert json.load(f)["seed"] == 7


def test_resolve_settings_precedence():
    s = interpret.resolve_settings(argparse.Namespace(shap_seed=7, seed=1, shap_precision="default",
                                                      precision="bf16", shap_batch_seqs=None, batch_seqs=3,
                                                      device="cpu"))
    assert (s["seed"], s["precision"], s["batch_seqs"], s["device"]) == (7, "default", 3, "cpu")
    with pytest.warns(UserWarning, match="bf16"):
        s = interpret.resolve_settings(argparse.Namespace(precision="bf16"))
    assert s["precision"] == "auto" and s["seed"] == 1234 and s["device"] == "auto"
    with pytest.raises(ValueError):
        interpret.resolve_settings(argparse.Namespace(shap_precision="bf16"))
    s = interpret.resolve_settings(argparse.Namespace())
    assert s == {"seed": 1234, "batch_seqs": None, "precision": "auto", "device": "auto", "num_shuffles": 20,
                 "profile_weighting": "chrombpnet", "save_references": False, "references": None}


def test_device_gpu_without_gpu_fails():
    import jax
    if jax.default_backend() == "gpu":
        pytest.skip("a GPU is present")
    with pytest.raises(RuntimeError, match="gpu"):
        interpret.device_context("gpu")


def test_cli_parser(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["interpret.py", "-g", "g.fa", "-r", "r.bed", "-m", "m.h5", "-o", "out",
                                      "-p", "profile", "--precision", "default", "--save-references"])
    a = interpret.fetch_interpret_args()
    assert a.profile_or_counts == ["profile"] and a.precision == "default" and a.save_references
    assert a.seed == 1234 and a.batch_seqs is None and a.device == "auto" and a.num_shuffles == 20


def test_importance_hdf5_to_bigwig(workdir, run_default):
    prefix = run_default.output_prefix
    out = str(workdir["dir"] / "counts")
    args = argparse.Namespace(hdf5=prefix + ".counts_scores.h5", regions=prefix + ".interpreted_regions.bed",
                              chrom_sizes=str(workdir["dir"] / "chrom.sizes"), output_prefix=out,
                              output_prefix_stats=out + ".stats.txt", debug_chr=None, tqdm=0)
    importance_hdf5_to_bigwig.main(args)
    proj = read(prefix + ".counts_scores.h5")["projected_shap"].sum(1)
    bw = pyBigWig.open(out + ".bw")
    centers = [cm for i, cm in enumerate(CENTERS) if i != EDGE]
    for i, (c, m) in enumerate(centers):
        np.testing.assert_array_equal(np.array(bw.values(c, m - L // 2, m + L // 2), np.float32),
                                      proj[i].astype(np.float32))
    bw.close()
    assert os.path.exists(out + ".stats.txt")


def test_legacy_style_file_readable(tmp_path):
    """A deepdish-like file (groups raw/shap/projected_shap, Blosc filter 32001) opens with h5py + hdf5plugin."""
    rng = np.random.RandomState(0)
    x = np.eye(4, dtype=np.int8)[rng.randint(0, 4, (3, 64))]
    hyp = rng.normal(size=(3, 64, 4))
    d = interpret.generate_shap_dict(x, hyp)
    path = str(tmp_path / "legacy.h5")
    with h5py.File(path, "w") as f:
        for key in scores_io.KEYS:
            f.create_group(key).create_dataset("seq", data=d[key]["seq"], **hdf5plugin.Blosc())
    assert h5py.h5z.filter_avail(32001)
    np.testing.assert_array_equal(scores_io.load_scores(path), d["projected_shap"]["seq"])
    np.testing.assert_array_equal(scores_io.load_scores(path, "shap"), hyp.transpose(0, 2, 1).astype(np.float16))


def test_import_has_no_side_effects():
    code = ("import sys\n"
            "for m in ('tensorflow', 'shap', 'deepdish', 'deeplift', 'tensorflow_probability'):\n"
            "    sys.modules[m] = None\n"
            "import chrombpnet.evaluation.interpret.interpret\n"
            "import chrombpnet.evaluation.make_bigwigs.importance_hdf5_to_bigwig\n"
            "assert 'keras' not in sys.modules and 'jax' not in sys.modules, 'heavy import at module level'\n"
            "import chrombpnet.evaluation.interpret.explainer\n")
    subprocess.run([sys.executable, "-c", code], check=True)
