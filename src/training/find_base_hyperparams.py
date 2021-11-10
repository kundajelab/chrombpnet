import pyfaidx
import math
import pandas as pd
import numpy as np
import pyBigWig
from modisco.visualization import viz_sequence
import matplotlib.pyplot as plt
import argparse
import splits
import os

def parse_args():
    parser=argparse.ArgumentParser(description="find hyper-parameters for chrombpnet")
    parser.add_argument("-b", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-o", "--out_dir", help="output dir for storing hyper-param TSV for bias and chrombpnet")
    parser.add_argument("-sr", "--negative-sampling-ratio", type=float, default=1.0, help="Ratio of negative to positive samples per epoch")
    parser.add_argument("-t", "--threshold-factor", type=float, default=0.5, help="A threshold is applied on maximum count of non-peak region for training bias model, which is set as this threshold x min(count over peak regions)")
    parser.add_argument("-fl", "--fold", type=str, required=True, help="Fold information - see splits.py to set folds")
    parser.add_argument("-ol", "--outputlen", type=int, required=True, help="Fold information - see splits.py to set folds")
    return parser.parse_args()



args = parse_args()

chroms_to_keep=splits.splits[int(args.fold)]["train"]
print(chroms_to_keep)

bw = pyBigWig.open(args.bigwig) 
peaks =  pd.read_csv(args.peaks,
                           sep='\t',
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
peaks = peaks[(peaks["chr"].isin(chroms_to_keep))]

nonpeaks =  pd.read_csv(args.nonpeaks,
                           sep='\t',
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
nonpeaks = nonpeaks[(nonpeaks["chr"].isin(chroms_to_keep))]



def get_cts(bw, peaks_df, width=1000):
    vals = []
    for i, r in peaks_df.iterrows():
        try:
            vals.append(np.nan_to_num(bw.values(r['chr'], 
                                                (r['start'] + r['summit']) - width//2,
                                                (r['start'] + r['summit']) + width//2)))
        except:
            #print(i)
            #print((r['start'] + r['summit'])- width//2, (r['start'] + r['summit']) + width//2)
            pass
    
    print(np.array(vals).shape)
    return np.sum(np.array(vals),axis=1)

peak_cnts = get_cts(bw, peaks, args.outputlen)
nonpeak_cnts = get_cts(bw, nonpeaks, args.outputlen)

print(peak_cnts.shape)

upper_thresh = np.quantile(nonpeak_cnts, 0.9999)
lower_thresh = np.quantile(nonpeak_cnts, 1-0.9999)
counts_threshold = np.min(peak_cnts)*args.threshold_factor
print(upper_thresh, counts_threshold)
upper_thresh = min(upper_thresh,counts_threshold)
counts_loss_weight = np.median(nonpeak_cnts[(nonpeak_cnts <= upper_thresh) & (nonpeak_cnts>=lower_thresh)])/10


file = open(os.path.join(args.out_dir, "bias_params.txt"),"w")
file.write("\t".join(["counts_sum_min_thresh", str(round(lower_thresh,2))]))
file.write("\n")
file.write("\t".join(["counts_sum_max_thresh", str(round(upper_thresh,2))]))
file.write("\n")
file.write("\t".join(["counts_loss_weight", str(round(counts_loss_weight,2))]))
file.write("\n")
file.write("\t".join(["filters", "128"]))
file.write("\n")
file.write("\t".join(["n_dil_layers", "4"]))
file.write("\n")
file.write("\t".join(["trainings_pts_post_thresh", str(sum((nonpeak_cnts <= upper_thresh) & (nonpeak_cnts>=lower_thresh)))]))
file.close()

if args.negative_sampling_ratio > 0:
    final_cnts = np.concatenate((peak_cnts,np.random.choice(nonpeak_cnts, replace=False, size=(int(args.negative_sampling_ratio*len(nonpeak_cnts))))))
else:
    final_cnts = peak_cnts

upper_thresh = np.quantile(final_cnts, 0.9999)
lower_thresh = np.quantile(final_cnts, 1-0.9999)
counts_loss_weight = np.median(nonpeak_cnts)/10


file = open(os.path.join(args.out_dir, "chrombpnet_params.txt"),"w")
file.write("\t".join(["cts_sum_min_thresh", str(round(lower_thresh,2))]))
file.write("\n")
file.write("\t".join(["cts_sum_max_thresh", str(round(upper_thresh,2))]))
file.write("\n")
file.write("\t".join(["counts_loss_weight", str(round(counts_loss_weight,2))]))
file.write("\n")
file.write("\t".join(["filters", "512"]))
file.write("\n")
file.write("\t".join(["n_dil_layers", "8"]))
file.write("\n")
file.write("\t".join(["trainings_pts_post_thresh", str(sum((final_cnts< upper_thresh) & (final_cnts>lower_thresh)))]))
file.close()
