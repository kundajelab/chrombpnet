"""chrombpnet.evaluation.modisco.run: modisco argv/env construction and a tiny real motifs + report run."""
import os
import subprocess
import sys

import h5py
import hdf5plugin
import numpy as np
import pytest

from chrombpnet.data import DefaultDataFile, get_default_data_path
from chrombpnet.evaluation.modisco import run
from chrombpnet.helpers.generate_reports import modisco_table


class Recorder:
    def __init__(self):
        self.calls = []

    def __call__(self, argv, **kwargs):
        self.calls.append((argv, kwargs))
        return subprocess.CompletedProcess(argv, 0)


@pytest.fixture
def recorder(monkeypatch):
    rec = Recorder()
    monkeypatch.setattr(run.subprocess, "run", rec)
    return rec


def fake_which(available):
    def which(name, mode=os.F_OK | os.X_OK, path=None):
        return "/opt/bin/" + name if name in available else None
    return which


def test_motifs_argv_reproduces_chrombpnet_1x_defaults(recorder, tmp_path):
    out = tmp_path / "sub" / "modisco.h5"
    assert run.modisco_motifs("scores.h5", str(out)) == str(out)
    (argv, kwargs), = recorder.calls
    assert os.path.basename(argv[0]) == "modisco"
    assert argv[1:] == ["motifs", "-i", "scores.h5", "-n", "50000", "-o", str(out), "-w", "500", "-l", "2",
                        "-z", "20", "-f", "5", "-t", "20", "-g", "5", "-j", "0"]
    assert kwargs["check"] is True
    assert (tmp_path / "sub").is_dir()


def test_motifs_argv_custom(recorder, tmp_path):
    run.modisco_motifs(tmp_path / "s.h5", tmp_path / "o.h5", max_seqlets=1000, window=400, n_leiden=3,
                       trim_size=30, initial_flank=10, threads=2)
    (argv, kwargs), = recorder.calls
    assert argv[1:] == ["motifs", "-i", str(tmp_path / "s.h5"), "-n", "1000", "-o", str(tmp_path / "o.h5"),
                        "-w", "400", "-l", "3", "-z", "20", "-f", "5", "-t", "30", "-g", "10", "-j", "0"]
    assert kwargs["env"]["NUMBA_NUM_THREADS"] == "2"
    assert kwargs["env"]["OMP_NUM_THREADS"] == "2"


def test_report_argv_with_meme_tomtom(recorder, monkeypatch, tmp_path):
    monkeypatch.setattr(run.shutil, "which", fake_which({"tomtom"}))
    meme = get_default_data_path(DefaultDataFile.motifs_meme)
    html = run.modisco_report("m.h5", str(tmp_path / "rep"), meme)
    (argv, kwargs), = recorder.calls
    assert argv[1:] == ["report-simple", "-i", "m.h5", "-o", str(tmp_path / "rep"), "-m", str(meme), "-n", "3"]
    assert kwargs["check"] is True
    assert html == os.path.join(str(tmp_path / "rep"), "motifs.html")
    assert (tmp_path / "rep").is_dir()


def test_report_argv_tomtom_lite_does_not_need_tomtom(recorder, monkeypatch, tmp_path):
    monkeypatch.setattr(run.shutil, "which", fake_which(set()))
    run.modisco_report("m.h5", str(tmp_path / "rep"), "db.meme", tomtom_lite=True)
    (argv, _), = recorder.calls
    assert argv[1:] == ["report-simple", "-i", "m.h5", "-o", str(tmp_path / "rep"), "-m", "db.meme", "-n", "3",
                        "-l"]


def test_report_without_tomtom_fails_early_and_suggests_lite(recorder, monkeypatch, tmp_path):
    monkeypatch.setattr(run.shutil, "which", fake_which(set()))
    with pytest.raises(RuntimeError, match="--tomtom-lite"):
        run.modisco_report("m.h5", str(tmp_path / "rep"), "db.meme")
    assert recorder.calls == []


def test_report_without_motif_database(recorder, monkeypatch, tmp_path):
    monkeypatch.setattr(run.shutil, "which", fake_which(set()))
    run.modisco_report("m.h5", str(tmp_path / "rep"), None)
    (argv, _), = recorder.calls
    assert argv[1:] == ["report-simple", "-i", "m.h5", "-o", str(tmp_path / "rep")]


def test_thread_env(monkeypatch, tmp_path):
    monkeypatch.delenv("NUMBA_NUM_THREADS", raising=False)
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.delenv("NUMBA_CACHE_DIR", raising=False)
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "6")
    env = run.modisco_env()
    assert env["NUMBA_NUM_THREADS"] == "6" and env["OMP_NUM_THREADS"] == "6"
    assert env["NUMBA_CACHE_DIR"] == str(tmp_path / "cache" / "chrombpnet" / "numba")
    assert os.access(env["NUMBA_CACHE_DIR"], os.W_OK)
    assert env["PATH"].split(os.pathsep)[0] == os.path.dirname(sys.executable) or \
        os.path.dirname(sys.executable) in env["PATH"].split(os.pathsep)

    assert run.modisco_env(threads=3)["NUMBA_NUM_THREADS"] == "3"

    monkeypatch.setenv("NUMBA_NUM_THREADS", "2")
    assert run.modisco_env()["NUMBA_NUM_THREADS"] == "2"
    assert run.modisco_env(threads=5)["NUMBA_NUM_THREADS"] == "5"

    monkeypatch.delenv("NUMBA_NUM_THREADS")
    monkeypatch.delenv("SLURM_CPUS_PER_TASK")
    assert int(run.modisco_env()["NUMBA_NUM_THREADS"]) == run.default_threads() >= 1

    monkeypatch.setenv("NUMBA_CACHE_DIR", str(tmp_path / "mine"))
    assert run.modisco_env()["NUMBA_CACHE_DIR"] == str(tmp_path / "mine")


def test_default_threads_honours_the_user_thread_settings(monkeypatch):
    # an explicit NUMBA_NUM_THREADS, else OMP_NUM_THREADS, is the budget; then the job's CPUs
    monkeypatch.delenv("NUMBA_NUM_THREADS", raising=False)
    monkeypatch.delenv("OMP_NUM_THREADS", raising=False)
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "6")
    assert run.default_threads() == 6

    monkeypatch.setenv("OMP_NUM_THREADS", "4")
    assert run.default_threads() == 4
    env = run.modisco_env()
    assert env["NUMBA_NUM_THREADS"] == env["OMP_NUM_THREADS"] == "4"

    monkeypatch.setenv("NUMBA_NUM_THREADS", "2")
    assert run.default_threads() == 2
    env = run.modisco_env()
    assert env["NUMBA_NUM_THREADS"] == env["OMP_NUM_THREADS"] == "2"
    env = run.modisco_env(threads=5)  # an explicit threads= still wins
    assert env["NUMBA_NUM_THREADS"] == env["OMP_NUM_THREADS"] == "5"

    monkeypatch.delenv("OMP_NUM_THREADS")
    for invalid in ("", "0", "-3", "two"):
        monkeypatch.setenv("NUMBA_NUM_THREADS", invalid)
        assert run.default_threads() == 6
        assert run.modisco_env()["NUMBA_NUM_THREADS"] == "6"
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "0")
    assert run.default_threads() >= 1


def test_report_threads(recorder, monkeypatch, tmp_path):
    monkeypatch.setattr(run.shutil, "which", fake_which(set()))
    monkeypatch.setenv("NUMBA_NUM_THREADS", "2")
    run.modisco_report("m.h5", str(tmp_path / "rep"), "db.meme", tomtom_lite=True, threads=3)
    run.modisco_report("m.h5", str(tmp_path / "rep"), "db.meme", tomtom_lite=True)
    (_, explicit), (_, default) = recorder.calls
    assert explicit["env"]["NUMBA_NUM_THREADS"] == explicit["env"]["OMP_NUM_THREADS"] == "3"
    assert default["env"]["NUMBA_NUM_THREADS"] == default["env"]["OMP_NUM_THREADS"] == "2"


def test_environment_is_not_mutated(recorder):
    before = dict(os.environ)
    run.modisco_motifs("s.h5", "o.h5", threads=1)
    assert dict(os.environ) == before


def test_missing_modisco_executable(monkeypatch):
    monkeypatch.setattr(run, "find_executable", lambda name: None)
    with pytest.raises(FileNotFoundError, match="modisco"):
        run.modisco_motifs("s.h5", "o.h5")


def test_modisco_failure_propagates(tmp_path):
    with pytest.raises(subprocess.CalledProcessError):
        run.modisco_motifs(str(tmp_path / "missing.h5"), str(tmp_path / "o.h5"), threads=1)


# ---------------------------------------------------------------- tiny real run

POSITIVE = ["AGATAAGA", "CCACCAGGGGGCGC"]   # GATA-like, CTCF-like
NEGATIVE = ["GGGACTTTCC"]                  # NF-kB-like, negative contributions


def write_scores(path, n=300, length=1000, seed=0):
    """Contribution-score h5 in the chrombpnet layout with planted motifs in the central 500 bp."""
    rng = np.random.RandomState(seed)
    idx = rng.randint(0, 4, size=(n, length))
    hyp = rng.normal(0, 0.02, size=(n, length, 4)).astype(np.float32)
    for i in range(n):
        for m, (motif, sign) in enumerate([(x, 1.0) for x in POSITIVE] + [(x, -1.0) for x in NEGATIVE]):
            if rng.rand() < 0.9:
                start = rng.randint(280 + 150 * m, 410 + 150 * m)
                for j, base in enumerate(motif):
                    b = "ACGT".index(base)
                    idx[i, start + j] = b
                    hyp[i, start + j, :] = -0.05 * sign
                    hyp[i, start + j, b] = sign
    raw = np.eye(4, dtype=np.int8)[idx].transpose(0, 2, 1)
    shap = hyp.transpose(0, 2, 1).astype(np.float16)
    with h5py.File(path, "w") as f:
        f.create_dataset("raw/seq", data=raw, **hdf5plugin.Blosc())
        f.create_dataset("shap/seq", data=shap, **hdf5plugin.Blosc())
        f.create_dataset("projected_shap/seq", data=(raw * shap).astype(np.float16), **hdf5plugin.Blosc())
    return path


def test_modisco_motifs_and_report_end_to_end(tmp_path):
    scores = write_scores(str(tmp_path / "chrombpnet_nobias.profile_scores.h5"))
    modisco_h5 = run.modisco_motifs(scores, str(tmp_path / "modisco_results_profile_scores.h5"),
                                    max_seqlets=2000, window=500, threads=min(4, run.default_threads()))
    with h5py.File(modisco_h5, "r") as f:
        assert f.attrs["window_size"] == 500
        assert len(f["pos_patterns"]) >= 2
        assert len(f["neg_patterns"]) >= 1
        n_pos = len(f["pos_patterns"])
        n_neg = len(f["neg_patterns"])
        width = f["pos_patterns"]["pattern_0"]["sequence"].shape[0]
    assert width == 20 + 2 * 5  # trim 20 + 2 x initial flank 5: the chrombpnet 1.x size (50 with modisco 2.5 defaults)

    report_dir = str(tmp_path / "evaluation" / "modisco_profile")
    html_path = run.modisco_report(modisco_h5, report_dir, get_default_data_path(DefaultDataFile.motifs_meme),
                                   tomtom_lite=True, threads=min(4, run.default_threads()))
    html = open(html_path).read()
    assert modisco_table.match_stat(html) == "pval"
    assert modisco_table.table_columns(html)[:4] == ["pattern", "num_seqlets", "modisco_cwm_fwd", "modisco_cwm_rev"]
    assert html.count("<td>pos_patterns.pattern_") == n_pos
    assert html.count("<td>neg_patterns.pattern_") == n_neg
    assert os.path.isfile(os.path.join(report_dir, "trimmed_logos", "pos_patterns.pattern_0.cwm.fwd.png"))
    upper = html.upper()
    assert "CTCF" in upper and "GATA" in upper

    kept = modisco_table.drop_negative_patterns(html)
    assert "neg_patterns" not in kept and kept.count("<td>pos_patterns.pattern_") == n_pos
    table, stat = modisco_table.load_motifs_table(html_path, "modisco_profile", drop_negative=True)
    assert stat == "pval"
    assert 'src="./modisco_profile/trimmed_logos/pos_patterns.pattern_0.cwm.fwd.png"' in table
    assert ">pos__0<" in table and "neg_" not in table  # chrombpnet 1.x report naming
