"""Contribution-score HDF5 files ({prefix}.{counts,profile}_scores.h5), the chrombpnet 1.x layout.

    /raw/seq             int8     (N, 4, L)  one-hot input
    /shap/seq            float16  (N, 4, L)  hypothetical contributions
    /projected_shap/seq  float16  (N, 4, L)  hypothetical contributions * one-hot input

chrombpnet 1.x wrote these with deepdish (PyTables Blosc, HDF5 filter 32001, blosclz level 9 + shuffle); this
writes the same filter with h5py + hdf5plugin, so `modisco motifs -i` (which imports hdf5plugin), finemo and
older readers keep working. Reading any of these files needs `import hdf5plugin` before opening them.
"""
import os

import h5py
import hdf5plugin
import numpy as np

KEYS = ("raw", "shap", "projected_shap")
DTYPES = {"raw": np.int8, "shap": np.float16, "projected_shap": np.float16}
CHUNK_ROWS = 16


def blosc_filter():
    return hdf5plugin.Blosc(cname="blosclz", clevel=9, shuffle=hdf5plugin.Blosc.SHUFFLE)


def score_arrays(seqs, hyp):
    """(n, L, 4) one-hot and hypothetical scores -> {key: (n, 4, L) array} as stored on disk."""
    seqs = np.asarray(seqs)
    hyp = np.asarray(hyp)
    assert seqs.shape == hyp.shape and seqs.shape[2] == 4, (seqs.shape, hyp.shape)
    return {
        "raw": np.transpose(seqs, (0, 2, 1)).astype(np.int8),
        "shap": np.transpose(hyp, (0, 2, 1)).astype(np.float16),
        "projected_shap": np.transpose(seqs * hyp, (0, 2, 1)).astype(np.float16),
    }


class ScoresWriter:
    """Stream (N, 4, L) score datasets into `path`; written to `path + ".partial"` and renamed on close()."""

    def __init__(self, path, n, seq_len):
        self.path = path
        self.tmp_path = path + ".partial"
        self.n = n
        self.written = 0
        self.f = h5py.File(self.tmp_path, "w", rdcc_nbytes=64 * 1024 * 1024)
        chunks = (max(1, min(n, CHUNK_ROWS)), 4, seq_len)
        self.ds = {}
        for key in KEYS:
            self.ds[key] = self.f.create_group(key).create_dataset(
                "seq", shape=(n, 4, seq_len), dtype=DTYPES[key], chunks=chunks, **blosc_filter())

    def write(self, start, seqs, hyp):
        arrays = score_arrays(seqs, hyp)
        stop = start + arrays["raw"].shape[0]
        for key in KEYS:
            self.ds[key][start:stop] = arrays[key]
        self.written += stop - start

    def close(self):
        self.f.close()
        if self.written != self.n:
            raise RuntimeError("{}: wrote {} of {} rows".format(self.tmp_path, self.written, self.n))
        os.replace(self.tmp_path, self.path)

    def abort(self):
        try:
            self.f.close()
        finally:
            if os.path.exists(self.tmp_path):
                os.remove(self.tmp_path)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc, tb):
        if exc_type is None:
            self.close()
        else:
            self.abort()


def save_scores(path, seqs, hyp):
    """Write a whole scores file at once."""
    seqs = np.asarray(seqs)
    with ScoresWriter(path, seqs.shape[0], seqs.shape[1]) as w:
        w.write(0, seqs, hyp)


def load_scores(path, key="projected_shap"):
    """Read one (N, 4, L) dataset of a scores file (new or chrombpnet 1.x deepdish/Blosc)."""
    with h5py.File(path, "r") as f:
        return f[key]["seq"][:]
