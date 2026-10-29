import hdf5plugin
import h5py
import argparse
from pathlib import Path
import pyfaidx

import numpy as np
import pandas as pd
from chrombpnet.training.utils.data_utils import one_hot

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


def compare_seqs(
    cur_chrom: str, cur_start: int, cur_stop: int, original_one_hot: np.array, genome
) -> None:
    """Get one-hot encoding for new chrom/start/stop, compare to
    one-hot encoding from original h5
    Args:
        cur_chrom: chromosome
        cur_start: start
        cur_stop: stop
        original_one_hot: one-hot written to h5
        genome: pyfaidx genome object
    """
    sequence = genome[cur_chrom][cur_start:cur_stop].seq
    new_one_hot = one_hot.dna_to_one_hot(sequence).transpose(2, 0, 1).squeeze()
    assert np.array_equal(new_one_hot, original_one_hot), (
        f"{cur_chrom}:{cur_start} - {cur_stop} doesn't match"
        f"original one hot encoding"
    )


def add_index_to_interpret_outputs(
    h5_in: Path,
    regions_in: Path,
    path_to_genome: str,
    inputlen: int = 2114,
    dry_run: bool = False,
    check_regions: bool = True,
) -> None:
    """Add chrom, start, stop indexes to h5 files. For updating chrombpnet v1.0 outputs

    Args:
        h5_in: path to h5
        regions_in: path to regions (chr, start, stop)
        path_to_genome: path to reference genome (make sure this is the same one used to
        make the original shap file)
        inputlen: for this version of chrombpnet, will always be 2114
        dry_run: don't write new fields, just check that sequences match interpret one-hots
        check_regions: only turn this to False if confident that using the correct regions file
    """
    with h5py.File(h5_in, mode="r") as f:
        # check to make sure that coords keys don't already exist
        keys = set(f.keys())
        if keys.intersection(['coords_chrom', 'coords_start', 'coords_end']):
            print('Keys have already been updated')
            return

        regions_df = pd.read_csv(
            regions_in,
            sep="\t",
            names=NARROWPEAK_SCHEMA,
            usecols=["chr", "start", "summit"],
        )
        # need to convert 1000 bp region back to 2114 region
        # use same procedure as in chrombpnet.evaluation.interpret.input_utils
        start = (regions_df["start"] + regions_df["summit"] - inputlen // 2).to_numpy()
        end = (regions_df["start"] + regions_df["summit"] + inputlen // 2).to_numpy()
        chroms = regions_df["chr"].to_list()

        # check that all one-hots match regions
        genome = pyfaidx.Fasta(path_to_genome, sequence_always_upper=True)
        if check_regions:
            one_hots = f["raw"]["seq"]
            for idx in range(0, len(chroms)):
                compare_seqs(
                    cur_chrom=chroms[idx],
                    cur_start=start[idx],
                    cur_stop=end[idx],
                    genome=genome,
                    original_one_hot=one_hots[idx, :, :],
                )

        if dry_run:
            print(f"Dry run complete. Regions {regions_in} matches h5 {h5_in}")
            return

    # write new values to original h5 file
    with h5py.File(h5_in, "a") as f:
        print(f.keys())
        # write new values to original h5 file
        f.create_dataset(name="coords_chrom", data=chroms)
        f.create_dataset(
            name="coords_start", data=start, dtype="int64", **hdf5plugin.Blosc()
        )
        f.create_dataset(
            name="coords_end", data=end, dtype="int64", **hdf5plugin.Blosc()
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Add chrom, start, stop to interpret files"
    )
    parser.add_argument("--h5_in", type=Path, help="Path to h5")
    parser.add_argument("--regions_in", type=Path, help="Path to regions bed file")
    parser.add_argument(
        "--path_to_genome", type=str, help="Path to genome reference fasta"
    )
    parser.add_argument(
        "--dry_run",
        type=bool,
        help="Run script without writing anything to file",
        default=False,
        required=False,
    )
    parser.add_argument(
        "--check_regions",
        type=bool,
        help="Check to make sure that regions matches h5",
        default=True,
        required=False,
    )

    args = parser.parse_args()

    add_index_to_interpret_outputs(
        h5_in=args.h5_in,
        regions_in=args.regions_in,
        path_to_genome=args.path_to_genome,
        dry_run=args.dry_run,
        check_regions=args.check_regions,
    )
