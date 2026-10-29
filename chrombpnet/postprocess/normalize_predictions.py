import pandas as pd
import numpy as np
import argparse
import pyBigWig
from pathlib import Path
from typing import Dict

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


def get_region_signal(bw, chrom: str, start: int, summit: int, window: int) -> np.array:
    """Utility function to get values from a given range in a bw
    Args:
        bw: pyBigWig object
        chrom: chromosome
        start: start
        summit: peak position
        window: number of bases before/after peak corresponding to predict region
    Returns:
        Bigwig values
    """

    values = bw.values(
        chrom, start + summit - window, start + summit + window, numpy=True
    )
    values = np.nan_to_num(values, nan=0.0)
    return values


def compute_normalization_constant(
    bw_file: Path, background_bed: Path, window: int = 500, pseudocount: int = 1
):
    """Calculate a dataset-wide constant based on the average signal in background regions
    Args:
        bw_file: path to bigwig
        background_bed: path to nonpeaks bed file
        window: how many bases plus/minus start does prediction window extend
        pseudocount: value added to counts to make log value valid

    Returns:
        Per-dataset value corresponding to average signal in background
    """
    with pyBigWig.open(str(bw_file)) as bw:
        background_regions = pd.read_csv(
            background_bed,
            sep="\t",
            header=None,
            names=NARROWPEAK_SCHEMA,
            usecols=["chr", "start", "summit"],
        )

        if len(background_regions) > 25000:
            background_regions = background_regions.sample(
                n=25000, replace=False, random_state=42
            )

        background_signals = np.empty(
            shape=background_regions.shape[0], dtype="float64"
        )
        counter = 0
        for row in background_regions.itertuples():
            background_values = get_region_signal(
                bw, row.chr, row.start, row.summit, window
            )
            background_signals[counter] = np.log(
                np.sum(background_values) + pseudocount
            )
            counter += 1

        background_mean = np.exp(np.mean(background_signals)) / (window * 2)
        return background_mean


def load_chrom_sizes(chrom_sizes_file):
    with open(chrom_sizes_file) as f:
        return [(line.split("\t")[0], int(line.strip().split("\t")[1])) for line in f]


def write_bigwig(
    in_bw_path,
    coords_chrom,
    coords_start,
    coords_end,
    chrom_sizes,
    output_path,
    normalization_constant,
):
    in_bw = pyBigWig.open(str(in_bw_path))
    bw = pyBigWig.open(str(output_path), "w")
    bw.addHeader(chrom_sizes)
    cur_chr = ""
    cur_end = 0

    chr_to_idx = {c: i for i, (c, _) in enumerate(chrom_sizes)}
    regions = [
        [
            coords_chrom[i],
            coords_start[i],
            coords_end[i],
            (coords_start[i] + coords_end[i]) // 2,
        ]
        for i in range(len(coords_chrom))
    ]
    # order_of_regs: indexes for order of regions
    order_of_regs = sorted(
        range(len(regions)),
        key=lambda x: (chr_to_idx.get(regions[x][0], 1e9), regions[x][1]),
    )
    for itr in range(len(order_of_regs)):
        i = order_of_regs[itr]
        i_chr, i_start, i_end, i_mid = regions[i]
        if i_chr != cur_chr:
            cur_chr = i_chr
            cur_end = 0

        if cur_end < i_start:
            cur_end = i_start
        assert i_end >= cur_end, f"Region end {i_end} < current end {cur_end}"

        next_end = i_end
        if itr + 1 < len(order_of_regs):
            n = order_of_regs[itr + 1]
            next_chr, next_start, _, next_mid = regions[n]
            if next_chr == i_chr and next_start < i_end:
                next_end = (i_mid + next_mid) // 2

        # note: unlike other instances where this fun is called, we are loading from a bigwig, not an h5
        # or a list of values/dataframe. So don't shift by i_start (want actual coordinate, not
        # an index relative to 0)

        # raw_vals = in_bw.values(i_chr, cur_end - i_start, next_end - i_start)
        raw_vals = in_bw.values(i_chr, cur_end, next_end)
        vals = np.array(raw_vals) / normalization_constant
        # vals = profile[i][cur_end - i_start:next_end - i_start]
        try:
            bw.addEntries(
                [i_chr] * len(vals),
                list(range(cur_end, next_end)),
                ends=list(range(cur_end + 1, next_end + 1)),
                values=vals.astype(float).tolist(),
            )
        except:
            print(f"Could not add {i_chr}, {cur_end}, {next_end} to bigwig")

        cur_end = next_end
    bw.close()
    in_bw.close()


def normalize_bigwigs(
    in_bw_obs: Path,
    in_bw_nobias: Path,
    in_bw_full: Path,
    in_bed: Path,
    out_dir: Path,
    in_nonpeaks: Path,
    chrom_sizes_in: Path,
) -> None:
    """Normalize bigwig with average signal in background regions. Get coordinates from associated h5
    Args:
        in_bw_obs/nobias/full: Path to bigwig to be normalized (observed, chrombpnet no bias pred, chrombpnet full pred)
        in_bed: Path to regions file for bw to be normalized
        out_dir: Where to save bigwigs
        in_nonpeaks: Path to bed file with nonpeak regions
        chrom_sizes_in: Path to chrom sizes file
    """
    normalization_constant = compute_normalization_constant(
        bw_file=in_bw_obs, background_bed=in_nonpeaks
    )
    print(normalization_constant)

    chrom_sizes = load_chrom_sizes(chrom_sizes_in)

    # dictionary with output files
    output_files: Dict = {
        in_bw_full: out_dir / "chrombpnet_pred_normalized.bw",
        in_bw_obs: out_dir / "obs_normalized.bw",
        in_bw_nobias: out_dir / "chrombpnet_nobias_pred_normalized.bw",
    }

    peaks = pd.read_csv(
        in_bed,
        sep="\t",
        header=None,
        names=NARROWPEAK_SCHEMA,
        usecols=["chr", "start", "summit"],
    )
    # use peaks for coordinates
    for in_bw in [in_bw_full, in_bw_nobias, in_bw_obs]:
        outpath = output_files[in_bw]
        write_bigwig(
            normalization_constant=normalization_constant,
            in_bw_path=in_bw,
            coords_chrom=peaks["chr"].to_numpy(),
            coords_start=(peaks["start"].to_numpy() + peaks["summit"].to_numpy() - 500),
            coords_end=(peaks["start"].to_numpy() + peaks["summit"].to_numpy() + 500),
            chrom_sizes=chrom_sizes,
            output_path=outpath,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Normalize bigwigs by dividing by average signal in background"
    )
    parser.add_argument("--in_bed", type=Path, required=True, help="peaks file")
    parser.add_argument(
        "--in_bw_obs", type=Path, required=True, help="observed bw to be normalized"
    )
    parser.add_argument(
        "--in_bw_nobias",
        type=Path,
        required=True,
        help="chrombpnet bias removed pred bw to be normalized",
    )
    parser.add_argument(
        "--in_bw_full",
        type=Path,
        required=True,
        help="chrombpnet full pred bw to be normalized",
    )
    parser.add_argument(
        "--out_dir", type=Path, required=True, help="Top dir to save normalized bw"
    )
    parser.add_argument(
        "--in_nonpeaks", type=Path, required=True, help="Bed file with nonpeaks"
    )
    parser.add_argument(
        "--chrom_sizes_in", type=Path, required=True, help="Path to chrom sizes file"
    )

    args = parser.parse_args()

    normalize_bigwigs(
        in_bw_obs=args.in_bw_obs,
        in_bw_nobias=args.in_bw_nobias,
        in_bw_full=args.in_bw_full,
        in_bed=args.in_bed,
        out_dir=args.out_dir,
        in_nonpeaks=args.in_nonpeaks,
        chrom_sizes_in=args.chrom_sizes_in,
    )
