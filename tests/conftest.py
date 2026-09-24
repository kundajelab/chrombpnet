"""Shared pytest fixtures.

* ``synthetic_data``: a tiny deterministic dataset (3 chromosomes, ~90 kb) with a FASTA + .fai, a chrom sizes
  file, 10-column narrowPeak peaks and nonpeaks, a per-base counts bigWig with signal at the peak summits and a
  fold JSON. Built once per session under a temporary directory; no network or reference data needed. Thin
  fixtures (``genome_fasta``, ``peaks_bed``, ...) return the individual paths as strings.
* ``goldens_dir`` / ``golden_path``: the legacy parity goldens produced by tests/goldens/make_goldens.py
  (pred_*.npz, shap_*_{counts,profile}.npz, fvals_*.npz, trace_*/...). Set CHROMBPNET_GOLDENS to their
  directory; tests using them are skipped otherwise.
* ``legacy_results`` / ``legacy_prefix``: a results tree with real models and outputs from chrombpnet 1.x, the
  one tests/goldens/run_goldens.sh reads. Tests using them are skipped when these variables are unset:

  - CHROMBPNET_LEGACY_RESULTS: the tree's root, laid out as
    ``bias_models/bias_model_<bias>/<prefix>/{models,evaluation,auxiliary}/`` (one directory per bias model) and
    ``full_models/<prefix>/{models,auxiliary}/``;
  - CHROMBPNET_LEGACY_PREFIX: the run name ``<prefix>``, which is both the run directory and the ``--file-prefix``
    of the bias-model files (``<prefix>_bias.h5``);
  - CHROMBPNET_LEGACY_BIAS: the bias model ``<bias>`` whose DeepSHAP goldens exist (the ``bias_<bias>`` in
    ``shap_bias_<bias>*.npz``); read by tests/parity/test_shap_parity.py and test_predict_parity.py.
"""
import os

# Before anything imports keras: a test that imports keras ahead of chrombpnet would otherwise get Keras'
# default TensorFlow backend.
os.environ.setdefault("KERAS_BACKEND", "jax")
os.environ.setdefault("MPLBACKEND", "Agg")

import json  # noqa: E402
from dataclasses import dataclass, field  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402
import pytest  # noqa: E402

SEED = 20260923
CHROM_SIZES = {"chr1": 32_000, "chr2": 30_000, "chr3": 28_000}
FOLD = {"test": ["chr3"], "valid": ["chr2"], "train": ["chr1"]}
INPUTLEN = 2114
MAX_JITTER = 500
# Summits stay this far from chromosome ends, so the hyperparameter scripts keep every region.
EDGE = INPUTLEN // 2 + MAX_JITTER
PEAK_SPACING = 2600
LEADING_N = 100          # telomere-like run of N at the start of each chromosome
SOFT_MASKED = ("chr2", 2_000, 2_500)   # lowercase stretch, as in UCSC/GENCODE FASTAs
GC_BLOCK = 1_000
MOTIF = "TGACTCA"        # AP-1 site planted at every peak summit


@dataclass
class SyntheticData:
    dir: Path
    genome: str
    fai: str
    chrom_sizes: str
    peaks: str
    nonpeaks: str
    bigwig: str
    fold: str
    sequences: dict = field(repr=False)      # chrom -> sequence (str) as written to the FASTA
    counts: dict = field(repr=False)         # chrom -> per-base counts (int32), as written to the bigWig
    peak_summits: list = field(repr=False)   # [(chrom, absolute summit)]
    nonpeak_summits: list = field(repr=False)


def make_sequences(rng, chrom_sizes=CHROM_SIZES):
    """Random sequence whose GC fraction changes every GC_BLOCK bases (0.25-0.75)."""
    seqs = {}
    for chrom, size in chrom_sizes.items():
        n_blocks = -(-size // GC_BLOCK)
        gc = np.repeat(rng.uniform(0.25, 0.75, n_blocks), GC_BLOCK)[:size]
        is_gc = rng.random(size) < gc
        strong = rng.random(size) < 0.5
        bases = np.where(is_gc, np.where(strong, "G", "C"), np.where(strong, "A", "T"))
        seqs[chrom] = bases
    return seqs


def place_summits(rng, chrom_sizes=CHROM_SIZES):
    """Peak summits every PEAK_SPACING bp (with a little jitter), and nonpeak centres half-way between them."""
    peaks, nonpeaks = [], []
    for chrom, size in chrom_sizes.items():
        centres = np.arange(EDGE + 200, size - EDGE - 200, PEAK_SPACING)
        summits = centres + rng.integers(-100, 101, len(centres))
        peaks += [(chrom, int(s)) for s in summits]
        mids = (summits[:-1] + summits[1:]) // 2
        nonpeaks += [(chrom, int(m)) for m in mids if EDGE <= m <= size - EDGE]
    return peaks, nonpeaks


def plant_peaks(rng, seqs, peak_summits):
    """GC-rich cores and an AP-1 motif at every peak summit; N at chromosome starts; one soft-masked stretch."""
    for chrom, summit in peak_summits:
        core = slice(summit - 150, summit + 150)
        seqs[chrom][core] = np.where(rng.random(300) < 0.65, rng.choice(["G", "C"], 300), rng.choice(["A", "T"], 300))
        start = summit - len(MOTIF) // 2
        seqs[chrom][start:start + len(MOTIF)] = list(MOTIF)
    out = {}
    for chrom, bases in seqs.items():
        bases = bases.copy()
        bases[:LEADING_N] = "N"
        if chrom == SOFT_MASKED[0]:
            lo, hi = SOFT_MASKED[1:]
            bases[lo:hi] = np.char.lower(bases[lo:hi])
        out[chrom] = "".join(bases.tolist())
    return out


def make_counts(rng, chrom_sizes, peak_summits):
    """Sparse background insertions plus a Gaussian-shaped pile-up of reads at every peak summit."""
    counts = {chrom: rng.poisson(0.02, size).astype(np.int32) for chrom, size in chrom_sizes.items()}
    for chrom, summit in peak_summits:
        width = int(rng.integers(60, 121))
        pos = np.arange(summit - 4 * width, summit + 4 * width)
        rate = rng.uniform(2.0, 10.0) * np.exp(-0.5 * ((pos - summit) / width) ** 2)
        counts[chrom][pos] += rng.poisson(rate).astype(np.int32)
    return counts


def narrowpeak_rows(rng, summits, prefix, peak_like):
    rows = []
    for i, (chrom, summit) in enumerate(summits):
        if peak_like:
            width = int(rng.integers(300, 701))
            offset = int(rng.integers(width // 4, 3 * width // 4))
            stats = rng.uniform(5, 50, 3).round(5)
            rows.append([chrom, summit - offset, summit - offset + width, "{}_{}".format(prefix, i), 1000, ".",
                         stats[0], stats[1], stats[2], offset])
        else:   # as written by `chrombpnet prep nonpeaks`: inputlen-wide, summit in the middle
            half = INPUTLEN // 2
            rows.append([chrom, summit - half, summit + half, ".", ".", ".", ".", ".", ".", half])
    return rows


def write_fasta(path, seqs, width=60):
    with open(path, "w") as f:
        for chrom, seq in seqs.items():
            f.write(">{}\n".format(chrom))
            for i in range(0, len(seq), width):
                f.write(seq[i:i + width] + "\n")
    import pyfaidx
    pyfaidx.Fasta(str(path)).close()   # writes the .fai next to the FASTA
    return str(path) + ".fai"


def write_bed(path, rows):
    with open(path, "w") as f:
        for row in rows:
            f.write("\t".join(str(v) for v in row) + "\n")


def write_bigwig(path, chrom_sizes, counts):
    """bedGraphToBigWig-like output: one interval per non-zero base, nothing stored for zeros."""
    import pyBigWig
    bw = pyBigWig.open(str(path), "w")
    bw.addHeader(list(chrom_sizes.items()))
    for chrom in chrom_sizes:
        idx = np.flatnonzero(counts[chrom])
        if len(idx):
            bw.addEntries([chrom] * len(idx), idx.tolist(), ends=(idx + 1).tolist(),
                          values=counts[chrom][idx].astype(float).tolist())
    bw.close()


def build_synthetic_dataset(directory, seed=SEED, chrom_sizes=CHROM_SIZES, fold=FOLD):
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    peak_summits, nonpeak_summits = place_summits(rng, chrom_sizes)
    sequences = plant_peaks(rng, make_sequences(rng, chrom_sizes), peak_summits)
    counts = make_counts(rng, chrom_sizes, peak_summits)

    genome = directory / "genome.fa"
    fai = write_fasta(genome, sequences)
    sizes = directory / "chrom.sizes"
    write_bed(sizes, [[c, s] for c, s in chrom_sizes.items()])
    peaks = directory / "peaks.narrowPeak"
    write_bed(peaks, narrowpeak_rows(rng, peak_summits, "peak", peak_like=True))
    nonpeaks = directory / "nonpeaks.bed"
    write_bed(nonpeaks, narrowpeak_rows(rng, nonpeak_summits, "nonpeak", peak_like=False))
    bigwig = directory / "counts.bw"
    write_bigwig(bigwig, chrom_sizes, counts)
    fold_path = directory / "fold_0.json"
    with open(fold_path, "w") as f:
        json.dump(fold, f, indent=4)

    return SyntheticData(dir=directory, genome=str(genome), fai=fai, chrom_sizes=str(sizes), peaks=str(peaks),
                         nonpeaks=str(nonpeaks), bigwig=str(bigwig), fold=str(fold_path), sequences=sequences,
                         counts=counts, peak_summits=peak_summits, nonpeak_summits=nonpeak_summits)


@pytest.fixture(scope="session")
def synthetic_data(tmp_path_factory):
    return build_synthetic_dataset(tmp_path_factory.mktemp("synthetic"))


@pytest.fixture(scope="session")
def genome_fasta(synthetic_data):
    return synthetic_data.genome


@pytest.fixture(scope="session")
def chrom_sizes_file(synthetic_data):
    return synthetic_data.chrom_sizes


@pytest.fixture(scope="session")
def peaks_bed(synthetic_data):
    return synthetic_data.peaks


@pytest.fixture(scope="session")
def nonpeaks_bed(synthetic_data):
    return synthetic_data.nonpeaks


@pytest.fixture(scope="session")
def counts_bigwig(synthetic_data):
    return synthetic_data.bigwig


@pytest.fixture(scope="session")
def fold_json(synthetic_data):
    return synthetic_data.fold


def _env_dir(var):
    value = os.environ.get(var)
    if not value:
        pytest.skip("{} is not set".format(var))
    path = Path(value)
    if not path.is_dir():
        pytest.skip("{}={} is not a directory".format(var, value))
    return path


@pytest.fixture(scope="session")
def goldens_dir():
    """Directory written by tests/goldens/run_goldens.sh (CHROMBPNET_GOLDENS)."""
    return _env_dir("CHROMBPNET_GOLDENS")


@pytest.fixture
def golden_path(goldens_dir):
    """golden_path("pred_chrombpnet.npz") -> Path inside goldens_dir; skips the test if that file is missing."""
    def get(name):
        path = goldens_dir / name
        if not path.exists():
            pytest.skip("golden {} not found in {}".format(name, goldens_dir))
        return path
    return get


@pytest.fixture(scope="session")
def legacy_results():
    """chrombpnet 1.x results tree with real bias / chrombpnet models (CHROMBPNET_LEGACY_RESULTS)."""
    return _env_dir("CHROMBPNET_LEGACY_RESULTS")


@pytest.fixture(scope="session")
def legacy_prefix():
    """Run name inside legacy_results (CHROMBPNET_LEGACY_PREFIX): full_models/<prefix>/, bias_models/*/<prefix>/."""
    value = os.environ.get("CHROMBPNET_LEGACY_PREFIX")
    if not value:
        pytest.skip("CHROMBPNET_LEGACY_PREFIX is not set")
    return value
