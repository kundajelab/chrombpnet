"""Generate parity goldens with the LEGACY chrombpnet stack (TF-Keras 2.x + kundajelab-shap).

Runs inside the legacy container (Python 3.8, TF 2.12, chrombpnet 1.x installed), CPU only. It imports the
container's own chrombpnet, never this checkout. Outputs float32/float64 arrays with gzip so the Keras 3 / JAX
port can be compared against them at tight tolerances.

Subcommands (each is a separate process because DeepSHAP needs TF1 graph mode):
    predict  logits/logcounts of legacy .h5 models on fixed sequences
    shap     DeepSHAP (kundajelab TFDeepExplainer) with injected, content-seeded dinucleotide-shuffled references
    fvals    model outputs on the DeepSHAP sequences and on each of their references
    trace    20 Adam steps from exported initial weights on fixed batches
"""
import argparse
import json
import os
import sys
import zlib

import numpy as np

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
INPUTLEN = 2114
OUTPUTLEN = 1000
NUM_SHUFS = 20


def log(msg):
    print(msg, flush=True)


def content_seed(onehot):
    return zlib.crc32(np.ascontiguousarray(onehot, dtype=np.int8).tobytes()) & 0xFFFFFFFF


def make_refs(seqs, num_shufs=NUM_SHUFS):
    from deeplift.dinuc_shuffle import dinuc_shuffle
    refs = np.empty((seqs.shape[0], num_shufs) + seqs.shape[1:], dtype=np.int8)
    for i, s in enumerate(seqs):
        refs[i] = dinuc_shuffle(s, num_shufs=num_shufs, rng=np.random.RandomState(content_seed(s)))
    return refs


def load_legacy_model(path):
    import tensorflow as tf
    from tensorflow.keras.models import load_model
    from tensorflow.keras.utils import get_custom_objects
    import chrombpnet.training.utils.losses as losses
    custom_objects = {"multinomial_nll": losses.multinomial_nll, "tf": tf}
    get_custom_objects().update(custom_objects)
    return load_model(path, compile=False)


def fetch_seqs(genome_path, chroms, centers):
    import pyfaidx
    from chrombpnet.training.utils.one_hot import dna_to_one_hot
    genome = pyfaidx.Fasta(genome_path)
    seqs = [str(genome[c][int(m) - INPUTLEN // 2:int(m) + INPUTLEN // 2]) for c, m in zip(chroms, centers)]
    genome.close()
    assert all(len(s) == INPUTLEN for s in seqs)
    return dna_to_one_hot(seqs).astype(np.int8)


def cmd_predict(a):
    import h5py
    os.makedirs(a.out, exist_ok=True)
    inputs_path = os.path.join(a.out, "pred_inputs.npz")
    if os.path.exists(inputs_path):
        d = np.load(inputs_path)
        seqs, chroms, centers, is_peak = d["seqs"], d["chrom"], d["center"], d["is_peak"]
    else:
        with h5py.File(a.coords_h5, "r") as f:
            chrom = np.array([c.decode() if isinstance(c, bytes) else c for c in f["coords/coords_chrom"][:]])
            center = f["coords/coords_center"][:]
            peak = f["coords/coords_peak"][:]
        rng = np.random.RandomState(1234)
        sel = []
        for flag in (1, 0):
            idx = np.where((chrom == a.chrom) & (peak == flag))[0]
            sel.append(rng.choice(idx, size=a.n // 2, replace=False))
        sel = np.sort(np.concatenate(sel))
        chroms, centers, is_peak = chrom[sel], center[sel], peak[sel]
        seqs = fetch_seqs(a.genome, chroms, centers)
        np.savez_compressed(inputs_path, seqs=seqs, chrom=chroms, center=centers, is_peak=is_peak)
    for spec in a.models:
        name, path = spec.split("=", 1)
        n = a.n_large if name in a.large else len(seqs)
        model = load_legacy_model(path)
        logits, logcounts = model.predict(seqs[:n].astype(np.float32), batch_size=32, verbose=0)
        np.savez_compressed(os.path.join(a.out, "pred_{}.npz".format(name)),
                            logits=logits.astype(np.float32), logcounts=logcounts.astype(np.float32), n=n)
        log("predict {}: {} rows, logcounts[:3]={}".format(name, n, logcounts[:3, 0]))


def read_regions(bed, n, chrom_groups, seed):
    import pandas as pd
    df = pd.read_csv(bed, sep="\t", names=NARROWPEAK_SCHEMA, dtype={"chr": str})
    rng = np.random.RandomState(seed)
    parts = []
    per = n // len(chrom_groups)
    for group in chrom_groups:
        sub = df[df["chr"].isin(group)]
        parts.append(sub.iloc[np.sort(rng.choice(len(sub), size=min(per, len(sub)), replace=False))])
    return pd.concat(parts).reset_index(drop=True)


def cmd_shap(a):
    import tensorflow as tf
    tf.compat.v1.disable_eager_execution()  # as legacy interpret.py does at import time
    import pyfaidx
    import shap
    import chrombpnet.evaluation.interpret.input_utils as input_utils
    import chrombpnet.evaluation.interpret.shap_utils as shap_utils

    os.makedirs(a.out, exist_ok=True)
    groups = [g.split(",") for g in a.chrom_groups]
    regions = read_regions(a.regions, a.n, groups, a.seed)
    genome = pyfaidx.Fasta(a.genome)
    seqs, used = input_utils.get_seq(regions, genome, INPUTLEN)
    genome.close()
    regions = regions[used].reset_index(drop=True)
    seqs = seqs.astype(np.int8)
    refs = make_refs(seqs)
    ref_by_key = {s.tobytes(): r for s, r in zip(seqs, refs)}

    model = load_legacy_model(a.model)
    np.savez_compressed(os.path.join(a.out, "shap_{}.npz".format(a.name)), x=seqs, refs=refs,
                        chrom=regions["chr"].values.astype(str),
                        start=regions["start"].values, summit=regions["summit"].values)

    for head in a.heads:
        captured = []

        def data(s, _refs=ref_by_key):
            return [_refs[np.ascontiguousarray(s[0], dtype=np.int8).tobytes()]]

        def combine(mult, orig_inp, bg_data, _cap=captured):
            _cap.append(np.array(mult[0], dtype=np.float32))
            return shap_utils.combine_mult_and_diffref(mult, orig_inp, bg_data)

        if head == "counts":
            target = tf.reduce_sum(model.outputs[1], axis=-1)
        else:
            target = shap_utils.get_weightedsum_meannormed_logits(model)
        explainer = shap.explainers.deep.TFDeepExplainer((model.input, target), data,
                                                         combine_mult_and_diffref=combine)
        log("shap {} {}: {} regions".format(a.name, head, len(seqs)))
        hyp = explainer.shap_values(seqs, progress_message=16)
        mult = np.stack(captured)
        assert mult.shape == (len(seqs), NUM_SHUFS, INPUTLEN, 4), mult.shape
        # one file per head keeps each file under 64 MB
        np.savez_compressed(os.path.join(a.out, "shap_{}_{}.npz".format(a.name, head)),
                            hyp=np.asarray(hyp, dtype=np.float64), mult=mult)


def cmd_fvals(a):
    d = np.load(os.path.join(a.out, "shap_{}.npz".format(a.name)))
    x, refs = d["x"], d["refs"]
    model = load_legacy_model(a.model)
    lx, cx = model.predict(x.astype(np.float32), batch_size=32, verbose=0)
    flat = refs.reshape((-1,) + refs.shape[2:]).astype(np.float32)
    lr, cr = model.predict(flat, batch_size=32, verbose=0)
    n, k = refs.shape[:2]
    np.savez_compressed(os.path.join(a.out, "fvals_{}.npz".format(a.name)),
                        logits_x=lx.astype(np.float32), logcounts_x=cx.astype(np.float32),
                        logits_refs=lr.reshape(n, k, -1).astype(np.float32),
                        logcounts_refs=cr.reshape(n, k, -1).astype(np.float32))
    log("fvals {}: done".format(a.name))


def cmd_trace(a):
    import pyBigWig
    import tensorflow as tf
    tf.config.experimental.enable_tensor_float_32_execution(False)
    import chrombpnet.training.models.bpnet_model as bpnet_model
    import chrombpnet.training.models.chrombpnet_with_bias_model as cbp_model

    d = np.load(os.path.join(a.out, "pred_inputs.npz"))
    n = a.steps * a.batch_size
    x = d["seqs"][:n]
    bw = pyBigWig.open(a.bigwig)
    y = np.stack([np.nan_to_num(np.array(bw.values(c, int(m) - OUTPUTLEN // 2, int(m) + OUTPUTLEN // 2)))
                  for c, m in zip(d["chrom"][:n], d["center"][:n])]).astype(np.float32)
    bw.close()
    logy = np.log(1 + y.sum(-1, keepdims=True)).astype(np.float32)

    class Args:
        seed = a.seed
        learning_rate = 1e-3

    params = {"filters": str(a.filters), "n_dil_layers": str(a.n_dil_layers),
              "counts_loss_weight": str(a.counts_loss_weight), "inputlen": str(INPUTLEN),
              "outputlen": str(OUTPUTLEN), "bias_model_path": a.bias_model or ""}
    arch = cbp_model if a.bias_model else bpnet_model
    model = arch.getModelGivenModelOptionsAndWeightInits(Args(), params)
    tdir = os.path.join(a.out, "trace_{}".format(a.name))
    os.makedirs(tdir, exist_ok=True)
    model.save(os.path.join(tdir, "init.h5"))
    steps = []
    for i in range(a.steps):
        sl = slice(i * a.batch_size, (i + 1) * a.batch_size)
        res = model.train_on_batch(x[sl].astype(np.float32), [y[sl], logy[sl]], return_dict=True)
        steps.append({k: float(v) for k, v in res.items()})
        log("trace {} step {}: {}".format(a.name, i, steps[-1]))
    model.save(os.path.join(tdir, "final.h5"))
    np.savez_compressed(os.path.join(tdir, "batches.npz"), x=x, y=y, logy=logy)
    with open(os.path.join(tdir, "trace.json"), "w") as f:
        json.dump({"steps": steps, "params": params, "batch_size": a.batch_size, "optimizer": "Adam(1e-3)",
                   "tf32": False}, f, indent=1)


def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd", required=True)
    for name in ("predict", "shap", "fvals", "trace"):
        s = sub.add_parser(name)
        s.add_argument("--out", required=True)
        s.add_argument("--genome")
    s = sub.choices["predict"]
    s.add_argument("--coords-h5", required=True)
    s.add_argument("--chrom", default="chr6")
    s.add_argument("--n", type=int, default=1000)
    s.add_argument("--n-large", type=int, default=256)
    s.add_argument("--large", nargs="*", default=["chrombpnet", "chrombpnet_nobias"])
    s.add_argument("--models", nargs="+", required=True, help="name=path.h5")
    for name in ("shap", "fvals"):
        s = sub.choices[name]
        s.add_argument("--name", required=True)
        s.add_argument("--model", required=True)
    s = sub.choices["shap"]
    s.add_argument("--regions", required=True)
    s.add_argument("--n", type=int, default=128)
    s.add_argument("--seed", type=int, default=1234)
    s.add_argument("--chrom-groups", nargs="+", default=["chr6", "chr21,chr22"])
    s.add_argument("--heads", nargs="+", default=["counts", "profile"])
    s = sub.choices["trace"]
    s.add_argument("--name", required=True)
    s.add_argument("--bigwig", required=True)
    s.add_argument("--filters", type=int, default=128)
    s.add_argument("--n-dil-layers", type=int, default=4)
    s.add_argument("--counts-loss-weight", type=float, default=10.0)
    s.add_argument("--bias-model", default=None)
    s.add_argument("--steps", type=int, default=20)
    s.add_argument("--batch-size", type=int, default=32)
    s.add_argument("--seed", type=int, default=1234)
    a = p.parse_args()
    {"predict": cmd_predict, "shap": cmd_shap, "fvals": cmd_fvals, "trace": cmd_trace}[a.cmd](a)


if __name__ == "__main__":
    sys.exit(main())
