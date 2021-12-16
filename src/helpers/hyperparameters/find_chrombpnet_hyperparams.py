import pyfaidx
import math
import pandas as pd
import numpy as np
import pyBigWig
from modisco.visualization import viz_sequence
import matplotlib.pyplot as plt
import argparse
from context import load_model_wrapper
from context import one_hot
import os
import json
import pyfaidx
import keras

def parse_args():
    parser=argparse.ArgumentParser(description="find hyper-parameters for chrombpnet")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-b", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
    parser.add_argument("-sr", "--negative-sampling-ratio", type=float, default=1.0, help="Ratio of negative to positive samples per epoch")
    parser.add_argument("-oth", "--outlier_threshold", type=float, default=0.9999, help="outlier threshold to use")
    parser.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")
    parser.add_argument("-il", "--inputlen", type=int, required=True, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, required=True, help="Prediction output length")
    parser.add_argument("-j", "--max_jitter", type=int, required=True, help="Sequence input length")
    parser.add_argument("-fil", "--filters", type=int, default=128, help="Number of filters to use in chrombpnet mode")
    parser.add_argument("-dil", "--n_dilation_layers", type=int, default=4, help="Number of dilation layers to use in chrombpnet model")
    parser.add_argument("-bmp", "--bias_model_path", type=str, required=True, help="path of bias model")
    parser.add_argument("-o", "--output_dir", help="output dir for storing hyper-param TSV for chrombpnet")
    return parser.parse_args()


def filter_edge_regions(peaks_df, bw, width, peaks_bool):
    """
    Filter regions in bed file that are on edges i.e regions that cannot be used to construct
    input length + jitter length of the sequence
    """
    input_shape = peaks_df.shape[0]

    # left edge case
    filtered = np.array((peaks_df['start'] + peaks_df['summit'] - width//2) < 0)
    peaks_df = peaks_df[~filtered]
    num_filtered = sum(filtered)

    # right edge case
    chrom_to_sizes = bw.chroms()
    filtered = []
    for i, r in peaks_df.iterrows():
        if r['start'] + r['summit'] + width//2 > chrom_to_sizes[r['chr']] :
            filtered.append(True)
        else:
            filtered.append(False)
    filtered=np.array(filtered)
    peaks_df = peaks_df[~filtered]
    num_filtered += sum(filtered)

    if peaks_bool:
        print("Number of peaks input: ",input_shape)
        print("Number of peaks filtered because the input/output is on the edge: ", num_filtered)
        print("Number of peaks being used: ",peaks_df.shape[0])
    else:
        print("Number of non peaks input: ",input_shape)
        print("Number of non peaks filtered because the input/output is on the edge: ", num_filtered)
        print("Number of non peaks being used: ",peaks_df.shape[0])

    return peaks_df

def get_seqs_cts(genome, bw, peaks_df, input_width=2114, output_width=1000):
    # output only counts - not log counts
    # output only one-hot encoded sequence
    vals = []
    seqs = []
    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - input_width//2):(r['start'] + r['summit'] + input_width//2)])
        seqs.append(sequence)
        bigwig_vals=np.nan_to_num(bw.values(r['chr'], 
                            (r['start'] + r['summit']) - output_width//2,
                            (r['start'] + r['summit']) + output_width//2))
        vals.append(bigwig_vals)
    return (np.sum(np.array(vals),axis=1), one_hot.dna_to_one_hot(seqs))

def adjust_bias_model_logcounts(bias_model, seqs, cts):
    """
    Given a bias model, sequences and associated counts, the function adds a 
    constant to the output of the bias_model's logcounts that minimises squared
    error between predicted logcounts and observed logcounts (infered from 
    cts). This simply reduces to adding the average difference between observed 
    and predicted to the "bias" (constant additive term) of the Dense layer.
    Typically the seqs and counts would correspond to training nonpeak regions.
    ASSUMES model_bias's last layer is a dense layer that outputs logcounts. 
    This would change if you change the model.
    """

    # safeguards to prevent misuse
    #assert(bias_model.layers[-1].name == "logcount_predictions")
    assert(bias_model.layers[-1].name == "logcounts")
    assert(bias_model.layers[-1].output_shape==(None,1))
    assert(isinstance(bias_model.layers[-1], keras.layers.Dense))

    print("Predicting within adjust counts")
    _, pred_logcts = bias_model.predict(seqs, verbose=True)
    delta = np.mean(np.log(1+cts) - pred_logcts.ravel())

    dw, db = bias_model.layers[-1].get_weights()
    bias_model.layers[-1].set_weights([dw, db+delta])
    return bias_model

if __name__=="__main__":

    args = parse_args()

    splits_dict=json.load(open(args.chr_fold_path))
    chroms_to_keep=splits_dict["train"]+splits_dict["valid"]
    test_chroms_to_keep=splits_dict["test"]
    print("evaluating hyperparameters on the following chromosomes",chroms_to_keep)

    bw = pyBigWig.open(args.bigwig) 
    genome = pyfaidx.Fasta(args.genome)

    # read peaks and non peaks    
    in_peaks =  pd.read_csv(args.peaks,
                           sep='\t',
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
    in_nonpeaks =  pd.read_csv(args.nonpeaks,
                           sep='\t',
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])

    # get train/valid peaks and test peaks seperately
    peaks = in_peaks[(in_peaks["chr"].isin(chroms_to_keep))]
    test_peaks = in_peaks[(in_peaks["chr"].isin(test_chroms_to_keep))]

    nonpeaks = in_nonpeaks[(in_nonpeaks["chr"].isin(chroms_to_keep))]
    test_nonpeaks = in_nonpeaks[(in_nonpeaks["chr"].isin(test_chroms_to_keep))]

    # step 1 filtering: filter peaks that are in the edges - which prevents us from making the inputlen regions
    peaks = filter_edge_regions(peaks, bw, args.inputlen+2*args.max_jitter, peaks_bool=1)
    test_peaks = filter_edge_regions(test_peaks, bw, args.inputlen, peaks_bool=1)
   
    nonpeaks = filter_edge_regions(nonpeaks, bw, args.inputlen+2*args.max_jitter, peaks_bool=0)
    test_nonpeaks = filter_edge_regions(test_nonpeaks, bw, args.inputlen, peaks_bool=0)

    # step 2 filtering: filter peaks that are outliers in train and valid set
    peak_cnts, peak_seqs = get_seqs_cts(genome, bw, peaks, args.inputlen, args.outputlen)
    nonpeak_cnts, nonpeak_seqs = get_seqs_cts(genome, bw, nonpeaks, args.inputlen, args.outputlen)

    if args.negative_sampling_ratio > 0:
        final_cnts = np.concatenate((peak_cnts,np.random.choice(nonpeak_cnts, replace=False, size=(int(args.negative_sampling_ratio*len(peak_cnts))))))
    else:
        final_cnts = peak_cnts

    upper_thresh = np.quantile(final_cnts, args.outlier_threshold)
    lower_thresh = np.quantile(final_cnts, 1-args.outlier_threshold)

    assert(len(peak_cnts) == peaks.shape[0])
    assert(len(nonpeak_cnts) == nonpeaks.shape[0])

    peaks = peaks[(peak_cnts< upper_thresh) & (peak_cnts>lower_thresh)]
    nonpeaks = nonpeaks[(nonpeak_cnts< upper_thresh) & (nonpeak_cnts>lower_thresh)]
    print("Number of peaks after removing outliers: ", peaks.shape[0])
    print("Number of nonpeaks after removing outliers: ", nonpeaks.shape[0])

    # store the filtered set of peaks
    # combine train and test peak set and store them in a new file
    frames = [peaks, test_peaks]
    all_peaks = pd.concat(frames)
    all_peaks.to_csv(os.path.join(args.output_dir, "filtered.peaks.bed"), sep="\t",  header=False, index=False)
    frames = [nonpeaks, test_nonpeaks]
    all_nonpeaks = pd.concat(frames)
    all_nonpeaks.to_csv(os.path.join(args.output_dir, "filtered.nonpeaks.bed"), sep="\t", header=False, index=False)

    # find counts loss weight for model training - using train and validation set
    counts_loss_weight = np.median(final_cnts[(final_cnts <= upper_thresh) & (final_cnts>=lower_thresh)])/10

    # adjust bias model for training  - using train and validation set
    bias_model = load_model_wrapper(args.bias_model_path)
    bias_model_scaled = adjust_bias_model_logcounts(bias_model, nonpeak_seqs[(nonpeak_cnts< upper_thresh) & (nonpeak_cnts>lower_thresh)], nonpeak_cnts[(nonpeak_cnts< upper_thresh) & (nonpeak_cnts>lower_thresh)])
    # save the new bias model
    bias_model_scaled.save(os.path.join(args.output_dir, "bias_model_scaled.h5"))

    # store the parameters being used  - using train and validation set
    file = open(os.path.join(args.output_dir, "chrombpnet_peak_params.txt"),"w")
    file.write("\t".join(["counts_sum_min_thresh", str(round(lower_thresh,2))]))
    file.write("\n")
    file.write("\t".join(["counts_sum_max_thresh", str(round(upper_thresh,2))]))
    file.write("\n")
    file.write("\t".join(["trainings_pts_post_thresh", str(sum((final_cnts< upper_thresh) & (final_cnts>lower_thresh)))]))
    file.write("\n")
    file.close()

    file = open(os.path.join(args.output_dir, "chrombpnet_model_params.txt"),"w")
    file.write("\t".join(["counts_loss_weight", str(round(counts_loss_weight,2))]))
    file.write("\n")
    file.write("\t".join(["filters", str(args.filters)]))
    file.write("\n")
    file.write("\t".join(["n_dil_layers", str(args.n_dilation_layers)]))
    file.write("\n")
    file.write("\t".join(["bias_model_path", os.path.join(args.output_dir, "bias_model_scaled.h5")]))
    file.write("\n")
    file.write("\t".join(["inputlen", str(args.inputlen)]))
    file.write("\n")
    file.write("\t".join(["max_jitter", str(args.max_jitter)]))
    file.write("\n")
    file.write("\t".join(["outputlen", str(args.outputlen)]))
    file.write("\n")
    file.write("\t".join(["chr_fold_path", str(args.chr_fold_path)]))
    file.close()
