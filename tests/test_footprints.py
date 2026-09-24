"""Marginal footprinting on tiny data (CPU) and the deepdish-compatible h5 layout of {prefix}_footprints.h5."""
import argparse
import json

import h5py
import numpy as np
import pandas as pd
import pytest

import chrombpnet  # noqa: F401  (sets KERAS_BACKEND before keras is imported)
import keras

from chrombpnet.evaluation.marginal_footprints import marginal_footprinting

INPUTLEN, OUTPUTLEN = 512, 256  # conv1 (21) + 2 dilated layers (2, 4) + profile conv (75): 512 -> 406 -> crop 256
CHROMS = [("chr1", 8000), ("chr2", 8000)]


def tiny_bpnet(filters=8, seed=0):
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
    return keras.Model(inputs=inp, outputs=[prof, counts])


def test_write_footprints_h5_layout(tmp_path):
    rng = np.random.RandomState(0)
    footprints = {}
    for motif in ["control", "tn5_1", "GATA+TAL"]:
        fp = rng.rand(OUTPUTLEN).astype(np.float32)
        footprints[motif] = [fp / fp.sum(), np.array([rng.rand() * 100], dtype=np.float32)]
    path = tmp_path / "x_footprints.h5"
    marginal_footprinting.write_footprints_h5(str(path), footprints)
    with h5py.File(path, "r") as f:
        assert f.attrs["DEEPDISH_IO_VERSION"] == 12
        assert sorted(f) == sorted(footprints)
        for motif, (fp, counts) in footprints.items():
            group = f[motif]
            assert isinstance(group, h5py.Group)
            assert group.attrs["TITLE"] == b"list:2"  # deepdish reads the group back as a 2-item list
            assert sorted(group) == ["i0", "i1"]
            assert group["i0"].dtype == np.float32 and group["i0"].shape == (OUTPUTLEN,)
            assert group["i1"].dtype == np.float32 and group["i1"].shape == (1,)
            np.testing.assert_array_equal(group["i0"][:], fp)
            np.testing.assert_array_equal(group["i1"][:], counts)


@pytest.fixture(scope="module")
def footprint_run(tmp_path_factory):
    d = tmp_path_factory.mktemp("footprints")
    rng = np.random.RandomState(0)
    with open(d / "genome.fa", "w") as f:
        for c, n in CHROMS:
            seq = "".join(rng.choice(list("ACGT"), n))
            f.write(">{}\n".format(c))
            f.writelines(seq[i:i + 60] + "\n" for i in range(0, n, 60))
    rows = [[c, m - 250, m + 250, ".", ".", ".", ".", ".", ".", 250] for c, _ in CHROMS for m in range(1000, 7001, 1000)]
    pd.DataFrame(rows).to_csv(d / "nonpeaks.bed", sep="\t", header=False, index=False)
    json.dump({"train": ["chr1"], "valid": [], "test": ["chr2"]}, open(d / "fold.json", "w"))
    pd.DataFrame([["tn5_1", "GCACAGTACAGAGCTG"], ["tn5_2", "GTGCACAGTTCTAGAGTGTGCAG"]]).to_csv(
        d / "motif_to_pwm.tsv", sep="\t", header=False, index=False)
    model = tiny_bpnet()
    model.save(str(d / "nobias.h5"))
    args = argparse.Namespace(model_h5=str(d / "nobias.h5"), regions=str(d / "nonpeaks.bed"), genome=str(d / "genome.fa"),
                              chr_fold_path=str(d / "fold.json"), output_prefix=str(d / "out" / "chrombpnet_nobias"),
                              motifs_to_pwm=str(d / "motif_to_pwm.tsv"), batch_size=4, ylim=[0.0, 0.05])
    (d / "out").mkdir()
    marginal_footprinting.main(args)
    return d, args, model


def test_footprints_outputs(footprint_run):
    d, args, model = footprint_run
    prefix = args.output_prefix
    for motif in ["control", "tn5_1", "tn5_2"]:
        assert (d / "out" / "chrombpnet_nobias.{}.footprint.png".format(motif)).exists()
    with h5py.File(prefix + "_footprints.h5", "r") as f:
        assert sorted(f) == ["control", "tn5_1", "tn5_2"]
        for motif in f:
            footprint, counts = f[motif]["i0"][:], f[motif]["i1"][:]
            assert footprint.shape == (OUTPUTLEN,) and counts.shape == (1,)
            np.testing.assert_allclose(footprint.sum(), 1.0, rtol=1e-5)  # mean of per-region normalized profiles
    response = open(prefix + "_max_bias_response.txt").read()
    kind, mean, values = response.split("_")  # parsed like this by make_html.py
    assert kind in ("corrected", "uncorrected")
    tn5 = [float(v) for v in values.split("/")]
    assert len(tn5) == 2 and float(mean) == round(np.mean(np.float32(tn5)), 3)
    assert kind == ("corrected" if max(tn5) < 0.003 else "uncorrected")


def test_footprint_matches_direct_prediction(footprint_run):
    d, args, model = footprint_run
    import pyfaidx
    from chrombpnet.training.utils.data_utils import get_seq
    regions = pd.read_csv(args.regions, sep="\t", names=marginal_footprinting.NARROWPEAK_SCHEMA)
    with pyfaidx.Fasta(args.genome) as genome:
        seqs = get_seq(regions[regions.chr == "chr2"], genome, INPUTLEN)  # test chromosomes only
    motif = "GCACAGTACAGAGCTG"
    x = seqs.copy()
    start = INPUTLEN // 2 - len(motif) // 2
    x[:, start:start + len(motif)] = np.eye(4, dtype=np.int8)[["ACGT".index(b) for b in motif]]
    fwd, rev = model.predict(x, verbose=0), model.predict(x[:, ::-1, ::-1], verbose=0)
    softmax = lambda z: np.exp(z - z.max(1, keepdims=True)) / np.exp(z - z.max(1, keepdims=True)).sum(1, keepdims=True)
    tot = softmax(fwd[0]) * (np.exp(fwd[1]) - 1) + (softmax(rev[0]) * (np.exp(rev[1]) - 1))[:, ::-1]
    expected = (tot / tot.sum(1, keepdims=True)).mean(0)
    with h5py.File(args.output_prefix + "_footprints.h5", "r") as f:
        np.testing.assert_allclose(f["tn5_1"]["i0"][:], expected, rtol=1e-4, atol=1e-7)
        np.testing.assert_allclose(f["tn5_1"]["i1"][:], (np.exp(fwd[1]) - 1 + np.exp(rev[1]) - 1).mean(0), rtol=1e-4)
