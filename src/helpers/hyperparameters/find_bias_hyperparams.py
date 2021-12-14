import pyfaidx
import math
import pandas as pd
import numpy as np
import pyBigWig
from modisco.visualization import viz_sequence
import matplotlib.pyplot as plt
import argparse
import os
import json

# input filters, dilation layers and outlier num
def parse_args():
    parser=argparse.ArgumentParser(description="find hyper-parameters for chrombpnet")
    parser.add_argument("-b", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-o", "--output_dir", help="output dir for storing hyper-param TSV for bias model")
    parser.add_argument("-t", "--bias_threshold_factor", type=float, default=0.5, help="A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x min(count over peak regions)")
    parser.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Fold information - see splits.py to set folds")
    parser.add_argument("-il", "--inputlen", type=int, required=True, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, required=True, help="Prediction output length")
    parser.add_argument("-oth", "--outlier_threshold", type=float, default=0.9999, help="outlier threshold to use")
    parser.add_argument("-fil", "--filters", type=int, default=128, help="Number of filters to use in bias mode")
    parser.add_argument("-dil", "--n_dilation_layers", type=int, default=4, help="Number of dilation layers to use in bias model")

    return parser.parse_args()

def get_cts(bw, peaks_df, width=1000):
    vals = []
    for i, r in peaks_df.iterrows():
        try:
            vals.append(np.nan_to_num(bw.values(r['chr'], 
                                            (r['start'] + r['summit']) - width//2,
                                            (r['start'] + r['summit']) + width//2)))
        except:
            pass
    return np.sum(np.array(vals),axis=1)


if __name__=="__main__":

    args = parse_args()

    splits_dict=json.load(open(args.chr_fold_path))
    chroms_to_keep=splits_dict["train"]
    print("evaluating hyperparameters on the following chromosomes",chroms_to_keep)

    bw = pyBigWig.open(args.bigwig) 
    peaks =  pd.read_csv(args.peaks,
                           sep='\t',
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
    peaks = peaks[(peaks["chr"].isin(chroms_to_keep))]

    nonpeaks =  pd.read_csv(args.nonpeaks,
                           sep='\t',
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
    nonpeaks = nonpeaks[(nonpeaks["chr"].isin(chroms_to_keep))]

    peak_cnts = get_cts(bw, peaks, args.outputlen)
    nonpeak_cnts = get_cts(bw, nonpeaks, args.outputlen)

    # find bias model hyper-parameters

    print("number of peak regions provided: ", peak_cnts.shape[0])
    print("number of non peak regions provided: ", nonpeak_cnts.shape[0])

    counts_threshold = np.min(peak_cnts)*args.bias_threshold_factor
    nonpeak_cnts = nonpeak_cnts[nonpeak_cnts < counts_threshold]

    upper_thresh = np.quantile(nonpeak_cnts, args.outlier_threshold)
    lower_thresh = np.quantile(nonpeak_cnts, 1-args.outlier_threshold)
    counts_loss_weight = np.median(nonpeak_cnts[(nonpeak_cnts < upper_thresh) & (nonpeak_cnts>lower_thresh)])/10


    file = open(os.path.join(args.output_dir, "bias_params.txt"),"w")
    file.write("\t".join(["counts_sum_min_thresh", str(round(lower_thresh,2))]))
    file.write("\n")
    file.write("\t".join(["counts_sum_max_thresh", str(round(upper_thresh,2))]))
    file.write("\n")
    file.write("\t".join(["counts_loss_weight", str(round(counts_loss_weight,2))]))
    file.write("\n")
    file.write("\t".join(["filters", str(args.filters)]))
    file.write("\n")
    file.write("\t".join(["n_dil_layers", str(args.n_dilation_layers)]))
    file.write("\n")
    file.write("\t".join(["trainings_pts_post_thresh", str(sum((nonpeak_cnts < upper_thresh) & (nonpeak_cnts>lower_thresh)))]))
    file.write("\n")
    file.write("\t".join(["inputlen", str(args.inputlen)]))
    file.write("\n")
    file.write("\t".join(["outputlen", str(args.outputlen)]))
    file.close()
