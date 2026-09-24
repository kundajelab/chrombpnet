"""CLI parser defaults, --help, and the pipeline wiring (interpret hand-off, modisco calls) without models."""
import argparse
import json
import os
import shutil
import subprocess
import sys
import types

import numpy as np
import pandas as pd
import pytest

import chrombpnet.parsers as parsers
import chrombpnet.pipelines as pipelines

TRAIN_REQUIRED = ["-g", "g.fa", "-c", "c.sizes", "-ibam", "x.bam", "-o", "out", "-d", "ATAC", "-p", "p.bed",
                  "-n", "n.bed", "-fl", "fold.json"]
QC_REQUIRED = ["-bw", "x.bw", "-g", "g.fa", "-c", "c.sizes", "-o", "out", "-d", "ATAC", "-p", "p.bed", "-n", "n.bed",
               "-fl", "fold.json"]
MINIMAL_ARGV = {
    "pipeline": ["pipeline"] + TRAIN_REQUIRED + ["-b", "bias.h5"],
    "train": ["train"] + TRAIN_REQUIRED + ["-b", "bias.h5"],
    "qc": ["qc"] + QC_REQUIRED + ["-cm", "cbp.h5", "-cmb", "nb.h5"],
    "bias pipeline": ["bias", "pipeline"] + TRAIN_REQUIRED + ["-b", "0.5"],
    "bias train": ["bias", "train"] + TRAIN_REQUIRED + ["-b", "0.5"],
    "bias qc": ["bias", "qc"] + QC_REQUIRED + ["-bm", "bias.h5"],
    "prep nonpeaks": ["prep", "nonpeaks", "-g", "g.fa", "-o", "out", "-p", "p.bed", "-c", "c.sizes", "-fl", "f.json"],
    "prep splits": ["prep", "splits", "-op", "out", "-c", "c.sizes", "-tcr", "chr1", "-vcr", "chr2"],
    "pred_bw": ["pred_bw", "-bm", "bias.h5", "-r", "r.bed", "-g", "g.fa", "-c", "c.sizes", "-op", "out"],
    "contribs_bw": ["contribs_bw", "-m", "m.h5", "-r", "r.bed", "-g", "g.fa", "-c", "c.sizes", "-op", "out"],
    "footprints": ["footprints", "-m", "m.h5", "-r", "r.bed", "-g", "g.fa", "-fl", "f.json", "-op", "out",
                   "-pwm_f", "m.tsv"],
}

LEGACY_TRACKABLES = ['logcount_predictions_loss', 'loss', 'logits_profile_predictions_loss',
                     'val_logcount_predictions_loss', 'val_loss', 'val_logits_profile_predictions_loss']
# chrombpnet 1.x defaults, per subcommand
LEGACY_TRAIN = dict(outlier_threshold=0.9999, ATAC_ref_path=None, DNASE_ref_path=None, num_samples=10000, inputlen=2114,
                    outputlen=1000, seed=1234, epochs=50, early_stop=5, learning_rate=0.001,
                    trackables=LEGACY_TRACKABLES, architecture_from_file=None, file_prefix=None, html_prefix="./",
                    bsort=False, tmpdir=None, no_st=False, batch_size=64)
LEGACY_DEFAULTS = {
    "pipeline": dict(LEGACY_TRAIN, negative_sampling_ratio=0.1, filters=512, n_dilation_layers=8, max_jitter=500),
    "train": dict(LEGACY_TRAIN, negative_sampling_ratio=0.1, filters=512, n_dilation_layers=8, max_jitter=500),
    "bias pipeline": dict(LEGACY_TRAIN, filters=128, n_dilation_layers=4, max_jitter=0),
    "bias train": dict(LEGACY_TRAIN, filters=128, n_dilation_layers=4, max_jitter=0),
    "qc": dict(file_prefix=None, batch_size=64, html_prefix="./"),
    "bias qc": dict(file_prefix=None, batch_size=64, html_prefix="./"),
    "prep nonpeaks": dict(inputlen=2114, stride=1000, neg_to_pos_ratio_train=2, blacklist_regions=None, seed=1234),
    "prep splits": dict(),
    "pred_bw": dict(chrombpnet_model=None, chrombpnet_model_nb=None, output_prefix_stats=None, batch_size=64, tqdm=1,
                    debug_chr=None, bigwig=None),
    "contribs_bw": dict(profile_or_counts=["counts", "profile"], output_prefix_stats=None, tqdm=1, debug_chr=None),
    "footprints": dict(batch_size=64, ylim=None),
}
# new flags: their defaults are the legacy behaviour (Adam, constant lr, no EMA, TF-like default precision,
# 30K interpret subsample with seed 1234, modisco -n 50000 -w 500 with MEME tomtom)
TRAINING_FLAGS = dict(optimizer="adam", muon_lr=None, ema=False, ema_momentum=0.999, lr_schedule="constant",
                      precision="default")
INTERPRET_FLAGS = dict(shap_seed=1234, shap_batch_seqs=None, shap_precision="auto")
MODISCO_FLAGS = dict(interpret_subsample=30000, modisco_max_seqlets=50000, modisco_window=500, tomtom_lite=False)
NEW_DEFAULTS = {
    "pipeline": dict(TRAINING_FLAGS, device="auto", **INTERPRET_FLAGS, **MODISCO_FLAGS),
    "train": dict(TRAINING_FLAGS, device="auto", **INTERPRET_FLAGS, **MODISCO_FLAGS),
    "bias pipeline": dict(TRAINING_FLAGS, device="auto", **INTERPRET_FLAGS, **MODISCO_FLAGS),
    "bias train": dict(TRAINING_FLAGS, device="auto", **INTERPRET_FLAGS, **MODISCO_FLAGS),
    "qc": dict(device="auto", **INTERPRET_FLAGS, **MODISCO_FLAGS),
    "bias qc": dict(device="auto", **INTERPRET_FLAGS, **MODISCO_FLAGS),
    "contribs_bw": dict(INTERPRET_FLAGS),
}


@pytest.mark.parametrize("command", sorted(MINIMAL_ARGV))
def test_parser_defaults_reproduce_legacy(command):
    args = vars(parsers.read_parser(MINIMAL_ARGV[command]))
    expected = dict(LEGACY_DEFAULTS[command], **NEW_DEFAULTS.get(command, {}))
    assert {k: args[k] for k in expected} == expected
    if command not in NEW_DEFAULTS:
        assert not set(args) & (set(TRAINING_FLAGS) | set(INTERPRET_FLAGS) | set(MODISCO_FLAGS) | {"device"})


HELP_ARGV = [[], ["pipeline"], ["train"], ["qc"], ["bias"], ["bias", "pipeline"], ["bias", "train"], ["bias", "qc"],
             ["prep"], ["prep", "nonpeaks"], ["prep", "splits"], ["pred_bw"], ["contribs_bw"], ["footprints"]]


@pytest.mark.parametrize("argv", HELP_ARGV, ids=lambda a: " ".join(a) or "top")
def test_help(argv, capsys):
    with pytest.raises(SystemExit) as exc:
        parsers.read_parser(argv + ["--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert out.startswith("usage:")
    if argv == ["pipeline"]:
        for flag in ["--optimizer", "--ema", "--precision", "--device", "--interpret-subsample", "--shap-seed",
                     "--modisco-max-seqlets", "--tomtom-lite"]:
            assert flag in out
        shap_batch_help = ("default: chosen automatically from the model size and available GPU memory; halved on "
                           "out-of-memory")
        assert "".join(shap_batch_help.split()) in "".join(out.split())  # argparse wraps (also at hyphens)


def test_help_subprocess_imports_no_deep_learning_framework():
    code = ("import sys\n"
            "sys.modules['tensorflow'] = None\n"
            "sys.argv = ['chrombpnet', 'bias', 'pipeline', '--help']\n"
            "import chrombpnet.CHROMBPNET as c\n"
            "try:\n"
            "    c.main()\n"
            "except SystemExit as e:\n"
            "    assert e.code == 0\n"
            "assert 'keras' not in sys.modules and 'jax' not in sys.modules\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "--modisco-window" in r.stdout


@pytest.mark.parametrize("command", ["snp_score", "modisco_motifs"])
def test_dead_commands_are_gone(command, capsys):
    with pytest.raises(SystemExit) as exc:
        parsers.read_parser([command])
    assert exc.value.code == 2


def test_ylim_is_two_floats():
    args = parsers.read_parser(MINIMAL_ARGV["footprints"] + ["--ylim", "0", "0.8"])
    assert args.ylim == [0.0, 0.8]
    with pytest.raises(SystemExit):
        parsers.read_parser(MINIMAL_ARGV["footprints"] + ["--ylim", "(0,0.8)"])


def test_new_flags_parse():
    args = parsers.read_parser(MINIMAL_ARGV["bias pipeline"] + [
        "--optimizer", "muon", "--muon-lr", "0.01", "--ema", "--ema-momentum", "0.99", "--lr-schedule", "cosine",
        "--precision", "bf16", "--device", "gpu", "--interpret-subsample", "100", "--shap-seed", "7",
        "--shap-batch-seqs", "4", "--shap-precision", "default", "--modisco-max-seqlets", "1000",
        "--modisco-window", "400", "--tomtom-lite"])
    assert (args.optimizer, args.muon_lr, args.ema, args.ema_momentum, args.lr_schedule) == ("muon", 0.01, True, 0.99,
                                                                                            "cosine")
    assert (args.precision, args.device) == ("bf16", "gpu")
    assert (args.interpret_subsample, args.shap_seed, args.shap_batch_seqs, args.shap_precision) == (100, 7, 4,
                                                                                                      "default")
    assert (args.modisco_max_seqlets, args.modisco_window, args.tomtom_lite) == (1000, 400, True)
    for bad in (["--optimizer", "sgd"], ["--precision", "fp16"], ["--device", "tpu"], ["--shap-precision", "bf16"]):
        with pytest.raises(SystemExit):
            parsers.read_parser(MINIMAL_ARGV["pipeline"] + bad)


def test_interpret_args_do_not_inherit_training_seed_or_precision():
    args = parsers.read_parser(MINIMAL_ARGV["pipeline"] + ["-s", "99", "--precision", "bf16", "--shap-seed", "5",
                                                           "--shap-batch-seqs", "8"])
    out = pipelines.interpret_args(argparse.Namespace(), args)
    assert (out.seed, out.precision, out.batch_seqs) == (5, "auto", 8)
    assert (args.seed, args.precision) == (99, "bf16")
    # namespaces that downstream pipelines build themselves have none of the new attributes
    out = pipelines.interpret_args(argparse.Namespace(seed=7), argparse.Namespace())
    assert (out.seed, out.precision, out.batch_seqs) == (1234, "auto", None)


# ---------------------------------------------------------------- pipeline wiring with fake steps

def _install(monkeypatch, name, **attrs):
    module = types.ModuleType(name)
    for key, value in attrs.items():
        setattr(module, key, value)
    monkeypatch.setitem(sys.modules, name, module)
    parent, _, leaf = name.rpartition(".")
    if parent in sys.modules:
        monkeypatch.setattr(sys.modules[parent], leaf, module, raising=False)
    return module


@pytest.fixture
def fake_steps(monkeypatch):
    calls = []
    import chrombpnet.evaluation.interpret  # noqa: F401  (parent packages of the fakes)
    import chrombpnet.evaluation.marginal_footprints  # noqa: F401
    import chrombpnet.evaluation.modisco  # noqa: F401
    import chrombpnet.helpers.generate_reports  # noqa: F401
    import chrombpnet.training  # noqa: F401
    import chrombpnet.training.utils.model_io as model_io

    class FakeModel:
        input_shape = (None, 2114, 4)
        output_shape = [(None, 1000), (None, 1)]

    monkeypatch.setattr(model_io, "load_model_wrapper", lambda model_h5: calls.append(("load", model_h5)) or FakeModel())

    def predict_main(args):
        calls.append(("predict", args.model_h5, args.peaks, args.nonpeaks, args.output_prefix, args.inputlen,
                      args.outputlen))

    def footprints_main(args):
        calls.append(("footprints", args.model_h5, args.regions, args.output_prefix))
        open(args.output_prefix + "_footprints.h5", "w").close()

    def interpret_main(args):
        calls.append(("interpret", args.model_h5, args.regions, args.output_prefix, tuple(args.profile_or_counts),
                      args.seed, args.precision, args.batch_seqs))

    def modisco_motifs(scores_h5, output_h5, max_seqlets=50000, window=500, n_leiden=2, trim_size=20,
                       initial_flank=5, threads=None):
        calls.append(("motifs", scores_h5, output_h5, max_seqlets, window, threads))

    def modisco_report(modisco_h5, output_dir, meme_file, tomtom_lite=False):
        calls.append(("report", modisco_h5, output_dir, meme_file, tomtom_lite))

    _install(monkeypatch, "chrombpnet.training.predict", main=predict_main)
    _install(monkeypatch, "chrombpnet.evaluation.marginal_footprints.marginal_footprinting", main=footprints_main)
    _install(monkeypatch, "chrombpnet.evaluation.interpret.interpret", main=interpret_main)
    _install(monkeypatch, "chrombpnet.evaluation.modisco.run", modisco_motifs=modisco_motifs,
             modisco_report=modisco_report)
    _install(monkeypatch, "chrombpnet.evaluation.modisco.convert_html_to_pdf",
             main=lambda html, pdf: calls.append(("pdf", html, pdf)))
    for name in ("make_html", "make_html_bias"):
        _install(monkeypatch, "chrombpnet.helpers.generate_reports." + name,
                 main=lambda args, name=name: calls.append((name, args.input_dir, args.command)))
    return calls


def _qc_namespace(tmp_path, n_peaks, **extra):
    out = tmp_path / "out"
    (out / "auxiliary").mkdir(parents=True)
    (out / "evaluation").mkdir()
    peaks = pd.DataFrame({0: ["chr1"] * n_peaks, 1: np.arange(n_peaks) * 1000, 2: np.arange(n_peaks) * 1000 + 500})
    for name in ("fp_filtered.bias_peaks.bed", "fp_filtered.peaks.bed"):
        peaks.to_csv(out / "auxiliary" / name, sep="\t", header=False, index=False)
    # like an argparse.Namespace a downstream pipeline builds for pipelines.bias_model_qc: none of the new flags
    return argparse.Namespace(bigwig="x.bw", bias_model="bias.h5", chrombpnet_model="cbp.h5",
                              chrombpnet_model_nb="nb.h5", genome="g.fa", chrom_sizes="c.sizes", output_dir=str(out),
                              data_type="ATAC", peaks="p.bed", nonpeaks="n.bed", chr_fold_path="f.json",
                              file_prefix="fp", batch_size=64, html_prefix="./", cmd_bias="qc", cmd="qc", **extra)


def test_bias_qc_wiring_defaults(tmp_path, fake_steps):
    args = _qc_namespace(tmp_path, n_peaks=12)
    pipelines.bias_model_qc(args)
    o = args.output_dir
    sub = os.path.join(o, "auxiliary/interpret_subsample/")
    meme = str(pipelines.get_default_data_path(pipelines.DefaultDataFile.motifs_meme))
    assert fake_steps == [
        ("load", "bias.h5"),
        ("predict", "bias.h5", os.path.join(o, "auxiliary/fp_filtered.bias_peaks.bed"),
         os.path.join(o, "auxiliary/fp_filtered.bias_nonpeaks.bed"), os.path.join(o, "evaluation/fp_bias"), 2114, 1000),
        ("interpret", "bias.h5", os.path.join(o, "auxiliary/fp_30K_subsample_peaks.bed"), sub + "fp_bias",
         ("counts", "profile"), 1234, "auto", None),
        # legacy order: profile motifs + report, then counts; `modisco motifs -n 50000 -w 500`, MEME tomtom
        ("motifs", sub + "fp_bias.profile_scores.h5", sub + "fp_modisco_results_profile_scores.h5", 50000, 500, None),
        ("report", sub + "fp_modisco_results_profile_scores.h5", os.path.join(o, "evaluation/modisco_profile/"), meme,
         False),
        ("motifs", sub + "fp_bias.counts_scores.h5", sub + "fp_modisco_results_counts_scores.h5", 50000, 500, None),
        ("report", sub + "fp_modisco_results_counts_scores.h5", os.path.join(o, "evaluation/modisco_counts/"), meme,
         False),
        ("pdf", os.path.join(o, "evaluation/modisco_counts/motifs.html"), os.path.join(o, "evaluation/fp_bias_counts.pdf")),
        ("pdf", os.path.join(o, "evaluation/modisco_profile/motifs.html"),
         os.path.join(o, "evaluation/fp_bias_profile.pdf")),
        ("make_html_bias", o, "qc"),
    ]
    assert len(pd.read_csv(os.path.join(o, "auxiliary/fp_30K_subsample_peaks.bed"), sep="\t", header=None)) == 12


def test_bias_qc_wiring_new_flags(tmp_path, fake_steps, monkeypatch):
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")  # many CPUs: the two heads still run one after the other
    args = _qc_namespace(tmp_path, n_peaks=12, interpret_subsample=5, shap_seed=3, shap_precision="default",
                         shap_batch_seqs=2, seed=99, precision="bf16", modisco_max_seqlets=1000, modisco_window=300,
                         tomtom_lite=True)
    pipelines.bias_model_qc(args)
    interpret = [c for c in fake_steps if c[0] == "interpret"]
    assert [c[5:] for c in interpret] == [(3, "default", 2)]
    modisco = [c for c in fake_steps if c[0] in ("motifs", "report")]
    assert [(c[0], c[1].rsplit("/", 1)[1]) for c in modisco] == [
        ("motifs", "fp_bias.profile_scores.h5"), ("report", "fp_modisco_results_profile_scores.h5"),
        ("motifs", "fp_bias.counts_scores.h5"), ("report", "fp_modisco_results_counts_scores.h5")]
    assert [c[3:] for c in modisco if c[0] == "motifs"] == [(1000, 300, None)] * 2
    assert [c[4] for c in modisco if c[0] == "report"] == [True, True]
    sub_peaks = pd.read_csv(os.path.join(args.output_dir, "auxiliary/fp_30K_subsample_peaks.bed"), sep="\t", header=None)
    peaks = pd.read_csv(os.path.join(args.output_dir, "auxiliary/fp_filtered.bias_peaks.bed"), sep="\t", header=None)
    pd.testing.assert_frame_equal(sub_peaks, peaks.sample(5, random_state=1234).reset_index(drop=True))
    # the pdf/html steps run after both modisco heads finished
    assert [c[0] for c in fake_steps][-3:] == ["pdf", "pdf", "make_html_bias"]


def test_chrombpnet_qc_wiring(tmp_path, fake_steps):
    args = _qc_namespace(tmp_path, n_peaks=3)
    pipelines.chrombpnet_qc(args)
    o = args.output_dir
    sub = os.path.join(o, "auxiliary/interpret_subsample/")
    kinds = [c[0] for c in fake_steps]
    assert kinds == ["load", "predict", "footprints", "interpret", "motifs", "report", "pdf", "make_html"]
    assert fake_steps[2] == ("footprints", "nb.h5", "n.bed", os.path.join(o, "evaluation/fp_chrombpnet_nobias"))
    assert os.path.exists(os.path.join(o, "auxiliary/fp_chrombpnet_nobias_footprints.h5"))
    assert fake_steps[3] == ("interpret", "nb.h5", os.path.join(o, "auxiliary/fp_30K_subsample_peaks.bed"),
                             sub + "fp_chrombpnet_nobias", ("profile",), 1234, "auto", None)
    assert fake_steps[4] == ("motifs", sub + "fp_chrombpnet_nobias.profile_scores.h5",
                             sub + "fp_modisco_results_profile_scores.h5", 50000, 500, None)
    assert fake_steps[5][1:3] == (sub + "fp_modisco_results_profile_scores.h5",
                                  os.path.join(o, "evaluation/modisco_profile/"))
    motifs = pd.read_csv(os.path.join(o, "auxiliary/motif_to_pwm.tsv"), sep="\t", header=None)
    assert list(motifs[0]) == ["tn5_1", "tn5_2", "tn5_3", "tn5_4", "tn5_5"]


def test_modisco_failure_propagates(tmp_path, fake_steps, monkeypatch):
    def failing_motifs(scores_h5, output_h5, **kwargs):
        raise subprocess.CalledProcessError(1, ["modisco", "motifs"])

    monkeypatch.setattr(sys.modules["chrombpnet.evaluation.modisco.run"], "modisco_motifs", failing_motifs)
    with pytest.raises(subprocess.CalledProcessError):
        pipelines.bias_model_qc(_qc_namespace(tmp_path, n_peaks=2))
    assert "pdf" not in [c[0] for c in fake_steps]


@pytest.mark.parametrize("numba_threads", [None, "2"])
def test_modisco_heads_run_one_after_the_other_with_the_job_threads(tmp_path, monkeypatch, numba_threads):
    """chrombpnet 1.x ran the profile and counts heads one after the other; each modisco process gets the user's
    NUMBA_NUM_THREADS if set, else the CPUs allocated to the job (run.default_threads), never the host CPU count."""
    import threading
    import chrombpnet.evaluation.modisco.run as run
    monkeypatch.setattr(os, "cpu_count", lambda: 64)  # a big host ...
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "4")  # ... and a 4-CPU job
    if numba_threads is None:
        monkeypatch.delenv("NUMBA_NUM_THREADS", raising=False)
    else:
        monkeypatch.setenv("NUMBA_NUM_THREADS", numba_threads)
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba"))
    calls = []

    def fake_run(argv, **kwargs):
        env = kwargs["env"]
        calls.append((argv[1], os.path.basename(argv[argv.index("-i") + 1]), env["NUMBA_NUM_THREADS"],
                      env["OMP_NUM_THREADS"], threading.current_thread() is threading.main_thread()))

    monkeypatch.setattr(run.subprocess, "run", fake_run)
    jobs = [(str(tmp_path / "{}_scores.h5".format(head)), str(tmp_path / "modisco_{}.h5".format(head)),
             str(tmp_path / "report_{}".format(head))) for head in ("profile", "counts")]
    pipelines.run_modisco_heads(argparse.Namespace(tomtom_lite=True), jobs, "motifs.meme")
    threads = numba_threads or "4"
    assert calls == [("motifs", "profile_scores.h5", threads, threads, True),
                     ("report-simple", "modisco_profile.h5", threads, threads, True),
                     ("motifs", "counts_scores.h5", threads, threads, True),
                     ("report-simple", "modisco_counts.h5", threads, threads, True)]
    assert not hasattr(pipelines, "available_cpus")


# ---------------------------------------------------------------- CHROMBPNET.main preflight checks

def _argv(command, out):
    argv = list(MINIMAL_ARGV[command])
    argv[argv.index("-o") + 1] = str(out)
    return ["chrombpnet"] + argv


@pytest.fixture
def no_tomtom(tmp_path, monkeypatch):
    """No `tomtom` next to the interpreter or on PATH (the linux pixi environments ship one)."""
    bindir = tmp_path / "bin"
    bindir.mkdir()
    monkeypatch.setattr(sys, "executable", str(bindir / "python"))
    monkeypatch.setenv("PATH", str(bindir))
    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "numba"))
    return bindir


@pytest.mark.parametrize("command", ["pipeline", "qc", "bias pipeline", "bias qc"])
def test_missing_tomtom_fails_before_any_output(tmp_path, monkeypatch, no_tomtom, command):
    import chrombpnet.CHROMBPNET as cli
    out = tmp_path / "out"
    monkeypatch.setattr(sys, "argv", _argv(command, out))
    with pytest.raises(RuntimeError, match="MEME's `tomtom` was not found on PATH.*--tomtom-lite"):
        cli.main()
    assert not out.exists()


@pytest.mark.parametrize("command,extra", [("pipeline", ["--tomtom-lite"]), ("bias qc", ["--tomtom-lite"]),
                                           ("train", []), ("bias train", [])])
def test_tomtom_not_needed(no_tomtom, command, extra):
    import chrombpnet.CHROMBPNET as cli
    cli.check_model_runtime(parsers.read_parser(MINIMAL_ARGV[command] + extra))


def test_tomtom_next_to_the_interpreter_is_found(no_tomtom):
    import chrombpnet.CHROMBPNET as cli
    tomtom = no_tomtom / "tomtom"
    tomtom.write_text("#!/bin/sh\n")
    tomtom.chmod(0o755)
    cli.check_model_runtime(parsers.read_parser(MINIMAL_ARGV["bias pipeline"]))


def test_device_check_only_for_gpu_and_cpu(monkeypatch):
    import chrombpnet.CHROMBPNET as cli
    import chrombpnet.training.runtime as runtime
    requested = []
    monkeypatch.setattr(runtime, "assert_gpu_if_requested", requested.append)
    monkeypatch.setenv("JAX_PLATFORMS", "")  # restored after the test
    for device in ("auto", "gpu", "cpu"):
        cli.check_model_runtime(parsers.read_parser(MINIMAL_ARGV["train"] + ["--device", device]))
    cli.check_model_runtime(parsers.read_parser(MINIMAL_ARGV["pred_bw"]))  # no --device
    assert requested == ["gpu", "cpu"]
    assert os.environ["JAX_PLATFORMS"] == "cpu"


@pytest.mark.parametrize("device", ["auto", "cpu"])
def test_device_check_initialises_no_other_jax_backend(device):
    """--device auto leaves JAX alone until the first model load (no CUDA context on every visible GPU during
    preprocessing); --device cpu restricts JAX to the CPU platform before any backend is created."""
    code = ("import argparse\n"
            "import chrombpnet.CHROMBPNET as cli\n"
            "cli.check_model_runtime(argparse.Namespace(cmd='train', device={!r}))\n"
            "import jax\n"
            "from jax._src import xla_bridge\n"
            "if {!r} == 'auto':\n"
            "    assert not xla_bridge.backends_are_initialized()\n"
            "else:\n"
            "    assert jax.config.jax_platforms == 'cpu', jax.config.jax_platforms\n"
            "    assert list(xla_bridge.backends()) == ['cpu'], list(xla_bridge.backends())\n").format(device, device)
    env = {k: v for k, v in os.environ.items() if k != "JAX_PLATFORMS"}
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, env=env)
    assert r.returncode == 0, r.stdout + r.stderr


# ---------------------------------------------------------------- prep commands through CHROMBPNET.main

def test_prep_splits(tmp_path, monkeypatch):
    import chrombpnet.CHROMBPNET as cli
    sizes = tmp_path / "chrom.sizes"
    sizes.write_text("chr1\t100\nchr2\t100\nchr3\t100\nchrX\t100\n")
    monkeypatch.setattr(sys, "argv", ["chrombpnet", "prep", "splits", "-op", str(tmp_path / "fold_0"), "-c", str(sizes),
                                      "-tcr", "chr1", "-vcr", "chr2", "chrX"])
    cli.main()
    assert json.load(open(tmp_path / "fold_0.json")) == {"test": ["chr1"], "valid": ["chr2", "chrX"],
                                                          "train": ["chr3"]}


def _write_fasta(path, chroms, rng):
    with open(path, "w") as f:
        for name, length in chroms:
            seq = "".join(rng.choice(list("ACGT"), length))
            f.write(">{}\n".format(name))
            f.writelines(seq[i:i + 60] + "\n" for i in range(0, length, 60))


def _prep_nonpeaks_inputs(d, rng):
    d.mkdir(parents=True)
    chroms = [("chr1", 40000), ("chr2", 40000), ("chr3", 40000)]
    _write_fasta(d / "genome.fa", chroms, rng)
    (d / "chrom.sizes").write_text("".join("{}\t{}\n".format(c, n) for c, n in chroms))
    json.dump({"train": ["chr1"], "valid": ["chr2"], "test": ["chr3"]}, open(d / "fold.json", "w"))
    rows = [[c, m - 250, m + 250, ".", ".", ".", ".", ".", ".", 250] for c, _ in chroms for m in (10000, 25000)]
    pd.DataFrame(rows).to_csv(d / "peaks.bed", sep="\t", header=False, index=False)
    (d / "blacklist.bed").write_text("chr1\t30000\t31000\nchr2\t5000\t5500\n")
    return rows


def _run_prep_nonpeaks(monkeypatch, d, prefix, blacklist=False):
    import chrombpnet.CHROMBPNET as cli
    argv = ["chrombpnet", "prep", "nonpeaks", "-g", str(d / "genome.fa"), "-o", prefix, "-p", str(d / "peaks.bed"),
            "-c", str(d / "chrom.sizes"), "-fl", str(d / "fold.json"), "-il", "1000", "-st", "500"]
    if blacklist:
        argv += ["-br", str(d / "blacklist.bed")]
    monkeypatch.setattr(sys, "argv", argv)
    cli.main()
    return pd.read_csv(prefix + "_negatives.bed", sep="\t", header=None)


@pytest.mark.needs_cli
def test_prep_nonpeaks(tmp_path, monkeypatch):
    if shutil.which("bedtools") is None:
        pytest.skip("bedtools not on PATH (run inside the pixi environment)")
    rows = _prep_nonpeaks_inputs(tmp_path / "in", np.random.RandomState(0))
    prefix = str(tmp_path / "out" / "fold_0")
    (tmp_path / "out").mkdir()
    negatives = _run_prep_nonpeaks(monkeypatch, tmp_path / "in", prefix)
    assert negatives.shape == (2 * 2 + 2 * 2 + 2 * 1, 10)  # 2 negatives per train/valid peak, 1 per test peak
    assert (negatives[9] == 500).all() and (negatives[2] - negatives[1] == 1000).all()
    for _, neg in negatives.iterrows():
        for c, s, e in [(r[0], r[1], r[2]) for r in rows]:
            assert not (neg[0] == c and neg[1] < e + 500 and neg[2] > s - 500)  # outside the slopped peaks


@pytest.mark.needs_cli
def test_prep_nonpeaks_paths_with_spaces(tmp_path, monkeypatch):
    """No shell: inputs and an output prefix with spaces (or shell syntax) give the same files as plain paths; the
    shell used to redirect `> /data/my run/...` to `/data/my`."""
    if shutil.which("bedtools") is None:
        pytest.skip("bedtools not on PATH (run inside the pixi environment)")
    rows = _prep_nonpeaks_inputs(tmp_path / "in", np.random.RandomState(0))
    _prep_nonpeaks_inputs(tmp_path / "in put $(x)", np.random.RandomState(0))
    (tmp_path / "out").mkdir()
    (tmp_path / "my run").mkdir()
    plain = str(tmp_path / "out" / "fold_0")
    spaced = str(tmp_path / "my run" / "fold 0")
    _run_prep_nonpeaks(monkeypatch, tmp_path / "in", plain, blacklist=True)
    negatives = _run_prep_nonpeaks(monkeypatch, tmp_path / "in put $(x)", spaced, blacklist=True)
    assert not (tmp_path / "my").exists()
    aux = ["blacklist_slop.bed", "candidates.bed", "exclude.bed", "exclude_unmerged.bed", "foreground.gc.bed",
           "genomewide_gc.bed", "negatives.bed", "peaks_slop.bed"]
    assert (sorted(os.listdir(plain + "_auxiliary")) == sorted(os.listdir(spaced + "_auxiliary"))
            == sorted(aux + ["negatives_compared_with_foreground.png"]))
    for name in aux:
        with open(os.path.join(plain + "_auxiliary", name)) as a, open(os.path.join(spaced + "_auxiliary", name)) as b:
            assert a.read() == b.read(), name
    with open(plain + "_negatives.bed") as a, open(spaced + "_negatives.bed") as b:
        assert a.read() == b.read()
    exclude = pd.read_csv(spaced + "_auxiliary/exclude.bed", sep="\t", header=None)
    assert len(exclude) == len(rows) + 2  # the slopped peaks and blacklist regions, none overlapping
    for _, neg in negatives.iterrows():
        assert not (neg[0] == "chr1" and neg[1] < 31500 and neg[2] > 29500)  # outside the slopped blacklist


@pytest.mark.needs_cli
def test_bedtools_sort_failure_is_not_hidden_by_merge(tmp_path, monkeypatch):
    """`bedtools sort | bedtools merge` without pipefail: a failed sort left merge an empty stdin (exit 0)."""
    import chrombpnet.CHROMBPNET as cli
    bedtools = shutil.which("bedtools")
    if bedtools is None:
        pytest.skip("bedtools not on PATH (run inside the pixi environment)")
    (tmp_path / "in.bed").write_text("chr1\t50\t80\nchr1\t10\t60\nchr2\t5\t9\n")
    cli.bedtools_sort_merge(str(tmp_path / "in.bed"), str(tmp_path / "merged.bed"))
    assert (tmp_path / "merged.bed").read_text() == "chr1\t10\t80\nchr2\t5\t9\n"
    bindir = tmp_path / "bin"
    bindir.mkdir()
    fake = bindir / "bedtools"
    fake.write_text('#!/bin/sh\nif [ "$1" = sort ]; then echo "sort: no space left" >&2; exit 1; fi\n'
                    'exec "{}" "$@"\n'.format(bedtools))
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(bindir) + os.pathsep + os.environ["PATH"])
    with pytest.raises(subprocess.CalledProcessError) as exc:
        cli.bedtools_sort_merge(str(tmp_path / "in.bed"), str(tmp_path / "merged.bed"))
    assert exc.value.cmd[:2] == ["bedtools", "sort"] and exc.value.returncode == 1


def test_prep_nonpeaks_fails_loudly_without_bedtools(tmp_path, monkeypatch):
    import chrombpnet.CHROMBPNET as cli
    rng = np.random.RandomState(0)
    _write_fasta(tmp_path / "genome.fa", [("chr1", 5000)], rng)
    (tmp_path / "chrom.sizes").write_text("chr1\t5000\n")
    json.dump({"train": ["chr1"], "valid": [], "test": []}, open(tmp_path / "fold.json", "w"))
    pd.DataFrame([["chr1", 2000, 2500, ".", ".", ".", ".", ".", ".", 250]]).to_csv(tmp_path / "peaks.bed", sep="\t",
                                                                                    header=False, index=False)
    monkeypatch.setenv("PATH", str(tmp_path / "empty_bin"))
    monkeypatch.setattr(sys, "argv", ["chrombpnet", "prep", "nonpeaks", "-g", str(tmp_path / "genome.fa"),
                                      "-o", str(tmp_path / "x"), "-p", str(tmp_path / "peaks.bed"),
                                      "-c", str(tmp_path / "chrom.sizes"), "-fl", str(tmp_path / "fold.json"),
                                      "-il", "1000"])
    with pytest.raises(FileNotFoundError, match="bedtools"):
        cli.main()
