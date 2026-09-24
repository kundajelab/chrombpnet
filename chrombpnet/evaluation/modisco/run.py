"""Run TF-MoDISco (the `modisco` command of the modisco>=2.5 / tfmodisco-lite package) on ChromBPNet scores.

`modisco_motifs` reads a contribution-score file written by `chrombpnet.evaluation.interpret`
(/raw/seq and /shap/seq, both (N, 4, L)) and writes the modisco results h5; `modisco_report` turns that into
`<output_dir>/motifs.html` (+ trimmed_logos/ and TOMTOM match logos) with `modisco report-simple`, the report
that the ChromBPNet HTML/PDF reports embed.
"""
import os
import shlex
import shutil
import subprocess
import sys
import tempfile

# modisco-lite 2.0.7, which chrombpnet 1.x ran, hard-coded these in its CLI (sliding window 20, seqlet flank 5,
# no final flank) and used trim 20 / initial flank 5. modisco 2.5 exposes them with different defaults
# (-t 30 -g 10), so every value is passed explicitly to reproduce the chrombpnet 1.x motifs.
SLIDING_WINDOW_SIZE = 20
SEQLET_FLANK_SIZE = 5
FINAL_FLANK_TO_ADD = 0
N_MATCHES = 3


def find_executable(name):
    """`name` from the environment chrombpnet runs in (next to the interpreter), else from PATH."""
    candidate = os.path.join(os.path.dirname(sys.executable), name)
    if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
        return candidate
    return shutil.which(name)


def _positive_int_env(name):
    try:
        value = int(os.environ[name])
    except (KeyError, ValueError):
        return None
    return value if value > 0 else None


def default_threads():
    """Thread count for MoDISco when none is given: an explicit user NUMBA_NUM_THREADS, else OMP_NUM_THREADS,
    else the CPUs allocated to this job (SLURM_CPUS_PER_TASK, else the CPU affinity mask, else os.cpu_count())."""
    for name in ("NUMBA_NUM_THREADS", "OMP_NUM_THREADS"):
        threads = _positive_int_env(name)
        if threads is not None:
            return threads
    threads = _positive_int_env("SLURM_CPUS_PER_TASK")
    if threads is not None:
        return threads
    if hasattr(os, "sched_getaffinity"):
        return max(1, len(os.sched_getaffinity(0)))
    return os.cpu_count() or 1


def _writable_cache_dir():
    cache_home = os.environ.get("XDG_CACHE_HOME") or os.path.join(os.path.expanduser("~"), ".cache")
    user = str(os.getuid()) if hasattr(os, "getuid") else "user"
    for path in (os.path.join(cache_home, "chrombpnet", "numba"),
                 os.path.join(tempfile.gettempdir(), "chrombpnet-numba-" + user)):
        try:
            os.makedirs(path, exist_ok=True)
        except OSError:
            continue
        if os.access(path, os.W_OK):
            return path
    return None


def modisco_env(threads=None):
    """Environment for the modisco subprocess: numba/OpenMP thread count, a writable numba cache, and the
    interpreter's bin directory on PATH (so `tomtom` from the same environment is found).

    threads=None uses default_threads(), so a user's NUMBA_NUM_THREADS / OMP_NUM_THREADS wins over the CPU count;
    an explicit threads overrides both variables."""
    env = dict(os.environ)
    if threads is None:
        threads = default_threads()
    threads = max(1, int(threads))
    env["NUMBA_NUM_THREADS"] = str(threads)
    env["OMP_NUM_THREADS"] = str(threads)
    if not env.get("NUMBA_CACHE_DIR"):
        cache_dir = _writable_cache_dir()
        if cache_dir is not None:
            env["NUMBA_CACHE_DIR"] = cache_dir
    bindir = os.path.dirname(sys.executable)
    path = env.get("PATH", "")
    if bindir not in path.split(os.pathsep):
        env["PATH"] = bindir + (os.pathsep + path if path else "")
    return env


def _modisco():
    exe = find_executable("modisco")
    if exe is None:
        raise FileNotFoundError("The `modisco` command (PyPI package modisco>=2.5.2) was not found next to {} "
                                "or on PATH.".format(sys.executable))
    return exe


def _run(argv, env):
    print("MoDISco: {} numba threads (NUMBA_NUM_THREADS), NUMBA_CACHE_DIR={}".format(
        env["NUMBA_NUM_THREADS"], env.get("NUMBA_CACHE_DIR")))
    print(" ".join(shlex.quote(a) for a in argv), flush=True)
    subprocess.run(argv, check=True, env=env)


def modisco_motifs(scores_h5, output_h5, max_seqlets=50000, window=500, n_leiden=2, trim_size=20,
                   initial_flank=5, threads=None):
    """`modisco motifs -i scores_h5` with chrombpnet's settings; returns output_h5.

    max_seqlets caps the seqlets per metacluster (runtime grows ~quadratically with it); window is the width
    around the region centre used for motif discovery (chrombpnet uses 500 bp to avoid AT-rich nucleosome
    flank motifs); threads (None: default_threads(), i.e. the user's NUMBA_NUM_THREADS / OMP_NUM_THREADS, else
    the CPUs allocated to the job) sets NUMBA_NUM_THREADS / OMP_NUM_THREADS of the modisco process.
    """
    argv = [_modisco(), "motifs",
            "-i", str(scores_h5),
            "-n", str(max_seqlets),
            "-o", str(output_h5),
            "-w", str(window),
            "-l", str(n_leiden),
            "-z", str(SLIDING_WINDOW_SIZE),
            "-f", str(SEQLET_FLANK_SIZE),
            "-t", str(trim_size),
            "-g", str(initial_flank),
            "-j", str(FINAL_FLANK_TO_ADD)]
    os.makedirs(os.path.dirname(os.path.abspath(str(output_h5))), exist_ok=True)
    _run(argv, modisco_env(threads))
    return output_h5


def modisco_report(modisco_h5, output_dir, meme_file, tomtom_lite=False, threads=None):
    """`modisco report-simple`: writes <output_dir>/motifs.html (returned) with the top TOMTOM matches of each
    pattern in meme_file (None skips motif matching).

    By default the matches come from MEME's `tomtom` binary (Pearson distance, q-value columns qval0..2,
    identical to chrombpnet 1.x). tomtom_lite=True uses memelite's TOMTOM-lite instead: no MEME install
    needed and much faster, but Euclidean distance and p-value columns (pval0..2).

    threads sets NUMBA_NUM_THREADS / OMP_NUM_THREADS of the report process as in modisco_motifs (None:
    default_threads()), so a caller that splits its CPUs between concurrent heads can cap the report too.
    """
    env = modisco_env(threads)
    argv = [_modisco(), "report-simple",
            "-i", str(modisco_h5),
            "-o", str(output_dir)]
    if meme_file is not None:
        argv += ["-m", str(meme_file), "-n", str(N_MATCHES)]
        if tomtom_lite:
            argv.append("-l")
        elif shutil.which("tomtom", path=env["PATH"]) is None:
            raise RuntimeError(
                "MEME's `tomtom` was not found on PATH; it is needed to match MoDISco motifs against {}. "
                "Install MEME (bioconda `meme`, included in chrombpnet's linux pixi environments) or use "
                "TOMTOM-lite instead (--tomtom-lite; reports p-values instead of q-values).".format(meme_file))
    os.makedirs(str(output_dir), exist_ok=True)
    _run(argv, env)
    return os.path.join(str(output_dir), "motifs.html")
