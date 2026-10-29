import sys
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.dirname(SCRIPT_DIR))

import h5py
import numpy as np
import argparse
from typing import List, Optional
from pathlib import Path
from chrombpnet_qc.bigwig_helper import write_bigwig, read_chrom_sizes, get_regions
import os
import hdf5plugin

os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"


def mean_predictions(
    predictions_list_str: List[str], output_path, input_len: int = 1000
):
    """Average across prediction h5 files from multiple folds.

    Args:
    predictions_list_str: paths to prediction files to be averaged together
    output_path: where to save averaged output
    input_len: number of bases in predict region
    """
    pred_logcounts: List = []
    logits: List = []
    chrom: List = []
    center: List = []
    for predictions_h5 in predictions_list_str:
        if not Path(predictions_h5).is_file():
            print(f"*********WARNING!************\n{predictions_h5} does not exist!")
            continue
        with h5py.File(predictions_h5, "r") as f:
            print(f"loading file {predictions_h5}")
            coord_chrom = f["coords"]["coords_chrom"][()]  # dim: num regions (1D)
            coord_center = f["coords"]["coords_center"][()]  # dim: num regions (1D)
            order = np.argsort(
                [
                    str(coord_chrom[pos]) + "_" + str(coord_center[pos])
                    for pos in range(f["coords"]["coords_chrom"][()].shape[0])
                ]
            )  # dim: num regions (1D)
            # pred_logcounts: predicted log(counts) -- i.e. from counts head
            # add 1D arrays to each list
            # flatten since some h5s store logcounts as (N, 1) and others as (N,)
            pred_logcounts.append(f["predictions"]["logcounts"][()].reshape(-1)[order])
            chrom.append(coord_chrom[order])
            center.append(coord_center[order])
            # pred_prob dim: # regions x 1000 -- these are probabilities from profile head
            pred_prob = f["predictions"]["profs"][()][order]
            # append the log of probabilities --> this is the logit
            logits.append(np.log(pred_prob))

    # check to make sure that all elements in the chromosome and center are the same across arrays
    compare_chroms = [np.array_equal(a, b) for a in chrom for b in chrom]
    assert all(compare_chroms), "Not all chromosomes are equal across folds!"
    compare_centers = [np.array_equal(a, b) for a in center for b in center]
    assert all(compare_centers), "Not all centers are equal across folds!"

    print("Averaging across files...")
    # average across log counts; dim for pred_logcounts_mean is # regions (1D) ...
    pred_logcounts_mean = np.nanmean(np.array(pred_logcounts), axis=0)
    # ...and exponentiate to get counts
    pred_counts_mean = np.exp(pred_logcounts_mean)

    # then take mean of profile logits, per position, across five folds...
    logit_mean = np.nanmean(np.array(logits), axis=0)
    # ...and softmax to get the probabilities back. dim is # regions x 1000
    pred_profile_mean = np.exp(logit_mean)/np.sum(np.exp(logit_mean), axis=1, keepdims=True)

    print("Finished averaging across files.")

    print("Writing outputs to file...")
    # write it all to file
    with h5py.File(output_path, "w") as f:
        f.create_dataset(
            name="coords/coords_chrom",
            data=chrom[0],
            dtype=h5py.string_dtype(encoding="ascii"),
            compression="gzip",
        )
        print("Created coors/coords_chrom dataset")
        f.create_dataset(
            name="coords/coords_center",
            data=center[0],
            dtype="int64",
            **hdf5plugin.Blosc(),
        )
        print("Created coords/coords_center dataset")
        f.create_dataset(
            name="coords/coords_start_dset",
            data=center[0] - input_len // 2,
            dtype="int64",
            **hdf5plugin.Blosc(),
        )
        print("Created coords/coords_start dataset")
        f.create_dataset(
            name="coords/coords_end_dset",
            data=center[0] + input_len // 2,
            dtype="int64",
            **hdf5plugin.Blosc(),
        )
        print("Created coords/coords_end dataset")
        f.create_dataset(
            name="predictions/logcounts",
            data=pred_logcounts_mean,
            dtype=float,
            **hdf5plugin.Blosc(),
        )
        print("Created predictions/logcounts dataset")
        f.create_dataset(
            name="predictions/profs",
            data=pred_profile_mean,
            dtype=float,
            **hdf5plugin.Blosc(),
        )
        print("Created predictions/profs dataset")


def make_bigwig_from_h5(h5_in: str, chrom_sizes_in: str, bw_out: str = None) -> None:
    """Loads data to pass to write_bigwig
    Args:
        h5_in: path to h5 with predictions
        regions_in: path to regions bed file (can be bed.gz). Should be 10 cols and the one that
        was used to create the h5
        chrom_sizes_in: path to chrom.sizes file. Should include all chroms in regions_in
        bw_out: where to save bigwig
    """
    regions = []
    # load profile to display in bigwig
    with h5py.File(h5_in, "r") as f:
        average_predicted_profile_probs = f["predictions"]["profs"][()]
        # reshape to (N, 1) so it broadcasts against profs (N, 1000)
        average_predicted_logcounts = f["predictions"]["logcounts"][()].reshape(-1, 1)
        average_predicted_profile = np.multiply(average_predicted_profile_probs, np.exp(average_predicted_logcounts) - 1)
        output_len = average_predicted_profile.shape[1] / 2
        n_regions = len(f["coords"]["coords_chrom"])
        for idx in range(0, n_regions):
            regions.append(
                [
                    f["coords"]["coords_chrom"][idx].decode(),
                    int(f["coords"]["coords_start_dset"][idx]),
                    int(f["coords"]["coords_end_dset"][idx]),
                    f["coords"]["coords_center"][idx],
                ]
            )

    # get chrom_sizes
    chrom_sizes = read_chrom_sizes(fname=chrom_sizes_in)

    if bw_out is None:
        bw_out = f"{h5_in.rstrip('h5')}bw"
    write_bigwig(
        data=average_predicted_profile,
        regions=regions,
        gs=chrom_sizes,
        bw_out=bw_out,
        debug_chr=None,
        use_tqdm=True,
        outstats_file=None,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Take a mean of prediction h5s and save the outputs as both hdf5 format"
    )
    parser.add_argument(
        "--prediction_h5s", nargs="+", type=str, required=False, help="prediction_h5s"
    )
    parser.add_argument(
        "-o", "--avg_h5", type=str, help="Output h5 file", required=True
    )
    parser.add_argument(
        "--average_preds",
        type=str,
        default="Y",
        help="Average pred h5s (Y/N)",
        choices=["Y", "N"],
    )
    parser.add_argument(
        "--make_bigwig",
        type=str,
        default="Y",
        help="Make prediction bigwig (Y/N)",
        choices=["Y", "N"],
    )
    parser.add_argument(
        "--chrom_sizes", type=str, required=True, help="Chrom sizes file"
    )
    parser.add_argument(
        "--bw_out", type=str, required=False, help="Name of output bw", default=None
    )
    args = parser.parse_args()

    if args.average_preds == "Y":
        print("Averaging bigwigs...")
        mean_predictions(args.prediction_h5s, args.avg_h5)
    else:
        print("Bypassing average_preds step...")

    if args.make_bigwig == "Y":
        print(args.avg_h5)
        make_bigwig_from_h5(
            h5_in=args.avg_h5, chrom_sizes_in=args.chrom_sizes, bw_out=args.bw_out
        )
    else:
        print("Bypassing make_bigwig step...")
