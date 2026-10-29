import argparse
import h5py
import pandas as pd
from pathlib import Path
from sklearn.metrics import roc_auc_score, average_precision_score


def calculate_binary_metrics(
    pred_in: Path, peaks_in: Path, neg_in: Path, output_dir: Path
) -> None:
    """Calculates binary metrics (AUPRC/AUROC) on predict files with test chromosomes
    Args:
        pred_in: path to prediction h5
        peaks_in: path to peaks file
        neg_in: path to negatives file
        output_dir: where to save file
     Returns:
         None
    """
    # load model
    f = h5py.File(pred_in, "r")

    # load peaks and non peaks
    peaks_df = pd.read_csv(peaks_in, header=None, sep="\t")
    # set value of fifth column to equal one (this will be the label column)
    peaks_df[4] = 1
    print(f"peaks_df.shape: {peaks_df.shape}")

    negs_df = pd.read_csv(neg_in, header=None, sep="\t")
    # set value of fifth column to equal zero (this will be the label column)
    negs_df[4] = 0
    print(f"negs_df.shape: {negs_df.shape}")

    # combine everything
    all_regions_df = pd.concat([peaks_df, negs_df]).reset_index(drop=True)
    all_regions_df.columns = [
        "chroms",
        "start",
        "end",
        "name",
        "label",
        "strand",
        "p",
        "q",
        "x",
        "summit",
    ]
    all_regions_df["coords_center"] = all_regions_df["start"] + all_regions_df["summit"]
    # peaks are shifted for some reason relative to prediction outputs
    all_regions_df.loc[all_regions_df["label"] == 1,"coords_center"] = all_regions_df.loc[all_regions_df["label"] == 1,"coords_center"] - 1057
    h5_df = pd.DataFrame(
        {
            "chroms": f["coords"]["coords_chrom"][()].astype("U8"),
            "coords_center": f["coords"]["coords_center"][()],
            "log_counts": f["predictions"]["logcounts"][()],
        }
    )
    # check that peaks mapped correctly
    coord_difference = h5_df['coords_center'][0] - all_regions_df['coords_center'][0]
    if coord_difference == 1057:
        print("Warning! Manually adjusting coords_center by 1057 for peaks regions.")
        all_regions_df.loc[all_regions_df['label'] == 1, 'coords_center'] = all_regions_df.loc[all_regions_df['label'] == 1, 'coords_center'] + 1057
    elif coord_difference == -1057:
        print("Warning! Manually adjusting coords_center by 1057 for peaks regions.")
        all_regions_df.loc[all_regions_df['label'] == 1, 'coords_center'] = all_regions_df.loc[all_regions_df['label'] == 1, 'coords_center'] - 1057
    print(sum(all_regions_df['coords_center'] - h5_df['coords_center']))
    h5_df.to_csv('h5_df.txt', sep=',', index=False)
    all_regions_df.to_csv('all_regions_df.5xt', sep=',', index=False)
    predictions = pd.merge(
        all_regions_df, h5_df, how="inner", on=["chroms", "coords_center"]
    )
    if predictions.shape[0] != all_regions_df.shape[0]:
        print("WARNING: predictions dataframe missing peaks and negatives regions.")

    auprc = average_precision_score(
        predictions["label"].astype(int), predictions["log_counts"]
    )
    auroc = roc_auc_score(predictions["label"].astype(int), predictions["log_counts"])
    auprc_baseline = sum(predictions["label"]) / len(predictions["label"])
    assert auprc_baseline == 0.2, "Incorrect value for auprc_baseline! Check that peaks regions match those in predictions."

    outfile = Path(f"{output_dir}/binary_metrics.txt")
    header = "auprc\tauroc\tauprc_baseline\n"
    out_str = f"{round(auprc, 3)}\t{round(auroc, 3)}\t{round(auprc_baseline, 3)}\n"
    with outfile.open("w") as f:
        f.write(f"{header}{out_str}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate AUPR/AUROC")
    parser.add_argument(
        "--h5_file", required=True, help="h5 file containing the predictions", type=Path
    )
    parser.add_argument("-p", "--peak_file", required=True, type=Path)
    parser.add_argument(
        "-n",
        "--neg_file",
        required=True,
        help="negative bed file; contains 5th column value 0",
        type=Path,
    )
    parser.add_argument(
        "-o", "--output_dir", required=True, help="directory for the output", type=Path
    )

    args = parser.parse_args()

    calculate_binary_metrics(
        pred_in=args.h5_file,
        peaks_in=args.peak_file,
        neg_in=args.neg_file,
        output_dir=args.output_dir,
    )
