import hdf5plugin
import h5py
import argparse
from pathlib import Path


def add_index_to_predict_outputs(
    h5_in: Path, inputlen: int = 1000
) -> None:
    """Add chrom, start, stop indexes to predict h5 files. For updating chrombpnet v1.0 outputs
    Args:
        h5_in: path to h5
        inputlen: for this version of chrombpnet, will always be 1000
    """
    # write new values to original h5 file
    with h5py.File(h5_in, "a") as f:
        coords_center = f['coords']['coords_center'][:]
        coords_start = coords_center - inputlen//2
        coords_end = coords_center + inputlen//2
        f.create_dataset(name="coords/coords_start_dset", data=coords_start, dtype="int64", **hdf5plugin.Blosc())
        f.create_dataset(name="coords/coords_end_dset", data=coords_end, dtype="int64", **hdf5plugin.Blosc())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add start, stop coords to predict files")
    parser.add_argument("--h5_in", type=Path, help="Path to h5")
    args = parser.parse_args()
    add_index_to_predict_outputs(h5_in=args.h5_in)
