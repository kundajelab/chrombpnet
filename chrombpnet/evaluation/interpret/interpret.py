# Adapted from chrombpnet-lite

import argparse
import contextlib
import json
import warnings

import numpy as np
import pandas as pd
import pyfaidx

import chrombpnet
import chrombpnet.evaluation.interpret.input_utils as input_utils
import chrombpnet.evaluation.interpret.scores_io as scores_io

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
PRECISIONS = ("auto", "highest", "high", "default")
DEVICES = ("auto", "gpu", "cpu")
PROFILE_WEIGHTINGS = ("chrombpnet", "chrombpnet_tf", "softmax_x", "tangermeme")

def fetch_interpret_args():
    parser = argparse.ArgumentParser(description="get sequence contribution scores for the model")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-m", "--model_h5", type=str, required=True, help="Path to trained model, can be both bias or chrombpnet model")
    parser.add_argument("-o", "--output-prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-d", "--debug_chr", nargs="+", type=str, default=None, help="Run for specific chromosomes only (e.g. chr1 chr2) for debugging")
    parser.add_argument("-p", "--profile_or_counts", nargs="+", type=str, default=["counts", "profile"], choices=["counts", "profile"],
                        help="use either counts or profile or both for running shap")
    parser.add_argument("--seed", type=int, default=1234, help="Seed for the dinucleotide-shuffled references (combined with each sequence's content, so scores do not depend on region order or batching)")
    parser.add_argument("--batch-seqs", type=int, default=None, help="Sequences per device step (each with its shuffled references); default: chosen from the model size and free device memory, halved on out-of-memory")
    parser.add_argument("--precision", type=str, default="auto", choices=PRECISIONS, help="Matmul/conv precision: 'auto' (default) = full float32 on CPU and TF32 on GPU (as chrombpnet 1.x); 'highest' = full float32 everywhere (very slow on some GPUs); 'default' = TF32 on Ampere+ GPUs")
    parser.add_argument("--device", type=str, default="auto", choices=DEVICES, help="'gpu' fails if JAX has no GPU; 'cpu' forces the CPU")
    parser.add_argument("--num-shuffles", type=int, default=20, help="Dinucleotide-shuffled references per sequence")
    parser.add_argument("--profile-weighting", type=str, default="chrombpnet", choices=PROFILE_WEIGHTINGS, help="Profile-head logit weights; 'chrombpnet' reproduces chrombpnet 1.x, the others are for comparison only")
    parser.add_argument("--save-references", action="store_true", help="Also save the references to {prefix}.references.npz (x, refs; large)")
    parser.add_argument("--references", type=str, default=None, help="npz with 'refs' (N, K, L, 4) (and optionally 'x') to use instead of generating references, e.g. from --save-references")

    args = parser.parse_args()
    return args


def generate_shap_dict(seqs, scores):
    assert(seqs.shape==scores.shape)
    assert(seqs.shape[2]==4)

    # construct a dictionary for the raw shap scores and the
    # the projected shap scores
    # MODISCO workflow expects one hot sequences with shape (None,4,inputlen)
    d = scores_io.score_arrays(seqs, scores)
    return {key: {'seq': d[key]} for key in scores_io.KEYS}


def _setting(args, name, default):
    # the chrombpnet CLI passes interpret options as shap_<name> next to same-named training options
    # (e.g. --seed, --precision), so the shap_ form wins
    for attr in ("shap_" + name, name):
        value = getattr(args, attr, None)
        if value is not None:
            return value, attr
    return default, None


def resolve_settings(args):
    seed, _ = _setting(args, "seed", 1234)
    batch_seqs, _ = _setting(args, "batch_seqs", None)
    precision, source = _setting(args, "precision", "auto")
    if precision not in PRECISIONS:
        if source == "precision":
            # a training-only value (e.g. bf16) reaching interpret through a shared namespace
            warnings.warn("DeepSHAP does not support precision {!r}; using 'auto'".format(precision))
            precision = "auto"
        else:
            raise ValueError("DeepSHAP precision must be one of {}, got {!r}".format(PRECISIONS, precision))
    device = getattr(args, "device", None) or "auto"
    if device not in DEVICES:
        raise ValueError("device must be one of {}, got {!r}".format(DEVICES, device))
    return {
        "seed": None if seed is None else int(seed),
        "batch_seqs": None if not batch_seqs else int(batch_seqs),
        "precision": precision,
        "device": device,
        "num_shuffles": int(getattr(args, "num_shuffles", None) or 20),
        "profile_weighting": getattr(args, "profile_weighting", None) or "chrombpnet",
        "save_references": bool(getattr(args, "save_references", False)),
        "references": getattr(args, "references", None),
    }


def device_context(device):
    import jax
    if device == "gpu" and jax.default_backend() != "gpu":
        raise RuntimeError("--device gpu requested but JAX runs on {!r} ({}); install the cuda13/cuda12 extra and "
                           "check the driver".format(jax.default_backend(), jax.devices()))
    if device == "cpu":
        return jax.default_device(jax.devices("cpu")[0])
    return contextlib.nullcontext()


def load_references(path, seqs):
    d = np.load(path)
    refs = d["refs"]
    if refs.ndim != 4 or refs.shape[0] != seqs.shape[0] or refs.shape[2:] != seqs.shape[1:]:
        raise ValueError("{}: refs have shape {}, expected ({}, K, {}, {})".format(
            path, refs.shape, seqs.shape[0], seqs.shape[1], seqs.shape[2]))
    if "x" in d.files and not np.array_equal(d["x"].astype(np.int8), seqs.astype(np.int8)):
        raise ValueError("{}: its sequences 'x' differ from the interpreted regions".format(path))
    return refs.astype(np.int8)


def interpret(model, seqs, output_prefix, profile_or_counts, seed=1234, batch_seqs=None, precision="auto",
              profile_weighting="chrombpnet", num_shuffles=20, save_references=False, references=None):
    from chrombpnet.evaluation.interpret.explainer import DeepLiftShap, timed_iter

    print("Seqs dimension : {}".format(seqs.shape))
    heads = [h for h in ("counts", "profile") if h in profile_or_counts]
    explainer = DeepLiftShap(model, heads=heads, profile_weighting=profile_weighting, precision=precision,
                             batch_seqs=batch_seqs, num_shuffles=num_shuffles)
    print("Generating {} shap scores".format(" and ".join("'{}'".format(h) for h in heads)))
    writers = {head: scores_io.ScoresWriter("{}.{}_scores.h5".format(output_prefix, head), seqs.shape[0],
                                            seqs.shape[1]) for head in heads}
    saved_refs = []
    additivity = {}
    try:
        chunks = explainer.iter_explain(seqs, references=references, seed=seed, additivity=additivity)
        for start, stop, refs, res in timed_iter(chunks, seqs.shape[0]):
            for head in heads:
                writers[head].write(start, seqs[start:stop], res[head]["hyp"])
            if save_references:
                saved_refs.append(refs)
    except BaseException:
        for w in writers.values():
            w.abort()
        raise
    for head in heads:
        # save the dictionary in HDF5 formnat
        print("Saving '{}' scores".format(head))
        writers[head].close()
        stat = additivity.get(head, {})
        print("'{}' summation-to-delta: max relative error {:.3g} over {} pairs".format(
            head, stat.get("max_rel_err", 0.0), stat.get("n_pairs", 0)))
    if save_references:
        refs = np.concatenate(saved_refs) if saved_refs else np.zeros((0, num_shuffles) + seqs.shape[1:], np.int8)
        np.savez_compressed("{}.references.npz".format(output_prefix), x=seqs.astype(np.int8), refs=refs)
    return additivity


def main(args):

    # check if the output directory exists
    #if not os.path.exists(os.path.dirname(args.output_prefix)):
    #    raise FileNotFoundError(errno.ENOENT, os.strerror(errno.ENOENT), os.path.dirname(args.output_prefix))

    import jax

    settings = resolve_settings(args)

    # write all the command line arguments to a json file, plus the DeepSHAP settings actually used
    record = dict(vars(args))
    record.update({k: v for k, v in settings.items() if k != "references" or v is not None})
    from chrombpnet.evaluation.interpret.explainer import resolve_precision
    record["precision_resolved"] = resolve_precision(settings["precision"])
    record.update({"jax_backend": jax.default_backend(), "jax_devices": [str(d) for d in jax.devices()],
                   "chrombpnet_version": chrombpnet.__version__})
    with open("{}.interpret.args.json".format(args.output_prefix), "w") as fp:
        json.dump(record, fp, ensure_ascii=False, indent=4, default=str)

    regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)

    if args.debug_chr:
        regions_df = regions_df[regions_df['chr'].isin(args.debug_chr)]

    with device_context(settings["device"]):
        model = input_utils.load_model_wrapper(args)
        # fail before fetching sequences or writing outputs if the model has no DeepSHAP rules here
        from chrombpnet.evaluation.interpret.explainer import check_model
        check_model(model, args.profile_or_counts)

        # infer input length
        inputlen = model.input_shape[1] # if bias model (1 input only)
        print("inferred model inputlen: ", inputlen)

        # load sequences
        # NOTE: it will pull out sequences of length inputlen
        #       centered at the summit (start + 10th column) and peaks used after filtering

        genome = pyfaidx.Fasta(args.genome)
        seqs, peaks_used = input_utils.get_seq(regions_df, genome, inputlen)
        genome.close()

        regions_df[peaks_used].to_csv("{}.interpreted_regions.bed".format(args.output_prefix), header=False, index=False, sep='\t')

        references = None
        if settings["references"]:
            references = load_references(settings["references"], seqs)

        interpret(model, seqs, args.output_prefix, args.profile_or_counts, seed=settings["seed"],
                  batch_seqs=settings["batch_seqs"], precision=settings["precision"],
                  profile_weighting=settings["profile_weighting"], num_shuffles=settings["num_shuffles"],
                  save_references=settings["save_references"], references=references)

if __name__ == '__main__':
    # parse the command line arguments
    args = fetch_interpret_args()
    main(args)
