import argparse
import pyfaidx
import pyBigWig
import pandas as pd
import numpy as np
import os
from chrombpnet.helpers.hyperparameters import param_utils as param_utils
from tensorflow import keras
import json

def parse_data_args():
    parser=argparse.ArgumentParser(description="find hyper-parameters for chrombpnet defined in src/training/models/chrombpnet_with_bias_model.py")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-i", "--bigwig", type=str, required=True, help="Bigwig of tn5 insertions. Ensure it is +4/-4 shifted")
    parser.add_argument("-p", "--peaks", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-n", "--nonpeaks", type=str, required=True, help="10 column bed file of non-peak regions, centered at summit (10th column)")
    parser.add_argument("-sr", "--negative_sampling_ratio", type=float, default=0.1, help="Ratio of negatives to positive samples per epoch")
    parser.add_argument("-oth", "--outlier-threshold", type=float, default=0.9999, help="threshold to use to filter outlies")
    parser.add_argument("-j", "--max-jitter", type=int, required=True, default=500, help="Maximum jitter applied on either side of region (default 500 for chrombpnet model)")
    parser.add_argument("-fl", "--chr-fold-path", type=str, required=True, help="Fold information - dictionary with test,valid and train keys and values with corresponding chromosomes")
    return parser

def parse_model_args(parser):
    # arguments here defined the following model - src/training/models/chrombpnet_with_bias_model.py
    parser.add_argument("-il", "--inputlen", type=int, required=True, help="Sequence input length")
    parser.add_argument("-ol", "--outputlen", type=int, required=True, help="Prediction output length")
    parser.add_argument("-fil", "--filters", type=int, default=512, help="Number of filters to use in chrombpnet mode")
    parser.add_argument("-dil", "--n-dilation-layers", type=int, default=8, help="Number of dilation layers to use in chrombpnet model")
    parser.add_argument("-b", "--bias-model-path", type=str, required=True, help="path of bias model")
    parser.add_argument("-op", "--output-prefix", help="output prefix for storing hyper-param TSV for chrombpnet")
    args = parser.parse_args()
    return args

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
    assert(bias_model.layers[-1].name == "logcounts" or bias_model.layers[-1].name == "logcount_predictions")
    assert(bias_model.layers[-1].output_shape==(None,1))
    assert(isinstance(bias_model.layers[-1], keras.layers.Dense))

    print("Predicting within adjust counts")
    _, pred_logcts = bias_model.predict(seqs, verbose=True)
    delta = np.mean(np.log(1+cts) - pred_logcts.ravel())

    dw, db = bias_model.layers[-1].get_weights()
    bias_model.layers[-1].set_weights([dw, db+delta])
    return bias_model


def main(args): 

    # read the fold information - we will evaluate hyperparams on the train+valid set and do nothing on the test set 
    splits_dict=json.load(open(args.chr_fold_path))
    chroms_to_keep=splits_dict["train"]+splits_dict["valid"]
    test_chroms_to_keep=splits_dict["test"]
    print("evaluating hyperparameters on the following chromosomes",chroms_to_keep)

    # read from bigwigw and fasta file
    bw = pyBigWig.open(args.bigwig) 
    genome = pyfaidx.Fasta(args.genome)

    # read peaks and non peaks    
    in_peaks =  pd.read_csv(args.peaks,
                           sep='\t',
                           header=None,
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
    in_nonpeaks =  pd.read_csv(args.nonpeaks,
                           sep='\t',
                           header=None,
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])

    assert(in_peaks.shape[0] != 0) # peaks file is empty
    assert(in_nonpeaks.shape[0] !=0) # non peaks file is empty
    assert(args.inputlen >= args.outputlen) # inputlen should be greater than the outputlen 
                                            # inputlen and outlen are chosen based on the filters and dilations layers used

    # get train/valid peaks and test peaks seperately
    peaks = in_peaks[(in_peaks["chr"].isin(chroms_to_keep))]
    test_peaks = in_peaks[(in_peaks["chr"].isin(test_chroms_to_keep))]

    nonpeaks = in_nonpeaks[(in_nonpeaks["chr"].isin(chroms_to_keep))]
    test_nonpeaks = in_nonpeaks[(in_nonpeaks["chr"].isin(test_chroms_to_keep))]

    # step 1 filtering: filter peaks that are in the edges - which prevents us from making the inputlen regions - do this for all splits
    peaks = param_utils.filter_edge_regions(peaks, bw, args.inputlen+2*args.max_jitter, peaks_bool=1)
    test_peaks = param_utils.filter_edge_regions(test_peaks, bw, args.inputlen, peaks_bool=1)
   
    nonpeaks = param_utils.filter_edge_regions(nonpeaks, bw, args.inputlen, peaks_bool=0)
    test_nonpeaks = param_utils.filter_edge_regions(test_nonpeaks, bw, args.inputlen, peaks_bool=0)

    # step 2 filtering: filter peaks that are outliers in train and valid set - no filtering on test set
    peak_cnts, _ = param_utils.get_seqs_cts(genome, bw, peaks, args.inputlen, args.outputlen)
    nonpeak_cnts, nonpeak_seqs = param_utils.get_seqs_cts(genome, bw, nonpeaks, args.inputlen, args.outputlen)    
    assert(len(peak_cnts) == peaks.shape[0])
    assert(len(nonpeak_cnts) == nonpeaks.shape[0])


    if args.negative_sampling_ratio > 0:
        final_cnts = np.concatenate((peak_cnts,np.random.choice(nonpeak_cnts, replace=False, size=(int(args.negative_sampling_ratio*len(peak_cnts))))))
    else:
        final_cnts = peak_cnts

    upper_thresh = np.quantile(final_cnts, args.outlier_threshold)
    lower_thresh = np.quantile(final_cnts, 1-args.outlier_threshold)

    peaks = peaks[(peak_cnts< upper_thresh) & (peak_cnts>lower_thresh)]
    nonpeaks = nonpeaks[(nonpeak_cnts< upper_thresh) & (nonpeak_cnts>lower_thresh)]

    print("Number of peaks after removing outliers: ", peaks.shape[0])
    print("Number of nonpeaks after removing outliers: ", nonpeaks.shape[0])

    # combine train valid and test peak set and store them in a new file
    if nonpeaks.shape[0] > peaks.shape[0]:
        train_nonpeaks = nonpeaks.sample(n=peaks.shape[0], random_state=1)
    else:
        train_nonpeaks = nonpeaks
        
    frames = [peaks, test_peaks]
    all_peaks = pd.concat(frames)
    all_peaks.to_csv("{}filtered.peaks.bed".format(args.output_prefix), sep="\t",  header=False, index=False)
    frames = [train_nonpeaks, test_nonpeaks]
    all_nonpeaks = pd.concat(frames)
    all_nonpeaks.to_csv("{}filtered.nonpeaks.bed".format(args.output_prefix), sep="\t", header=False, index=False)

    # find counts loss weight for model training - using train and validation set
    counts_loss_weight = np.median(final_cnts[(final_cnts <= upper_thresh) & (final_cnts>=lower_thresh)])/10
    assert(counts_loss_weight != 0)
    #assert(counts_loss_weight > 1.0) # counts loss weight can go less than 1.0 if you have very low-read depth - make sure you have enough density in counts
                                     # check peak-calling

    if counts_loss_weight < 1.0:
        counts_loss_weight = 1.0
        print("WARNING: you are training on low-read depth data")

    # adjust bias model for training  - using train and validation set
    # the bias model might be trained on a difference read depth compared to the given data - so this step scales the bias model to account for that
    bias_model = param_utils.load_model_wrapper(args.bias_model_path)
    bias_model_scaled = adjust_bias_model_logcounts(bias_model, nonpeak_seqs[(nonpeak_cnts< upper_thresh) & (nonpeak_cnts>lower_thresh)], nonpeak_cnts[(nonpeak_cnts< upper_thresh) & (nonpeak_cnts>lower_thresh)])
    # save the new bias model
    bias_model_scaled.save("{}bias_model_scaled.h5".format(args.output_prefix))

    # store the parameters being used  - in a TSV file
    file = open("{}chrombpnet_data_params.tsv".format(args.output_prefix),"w")
    file.write("\t".join(["counts_sum_min_thresh", str(round(lower_thresh,2))]))
    file.write("\n")
    file.write("\t".join(["counts_sum_max_thresh", str(round(upper_thresh,2))]))
    file.write("\n")
    file.write("\t".join(["trainings_pts_post_thresh", str(sum((final_cnts<upper_thresh) & (final_cnts>lower_thresh)))]))
    file.write("\n")
    file.close()

    file = open("{}chrombpnet_model_params.tsv".format(args.output_prefix),"w")
    file.write("\t".join(["counts_loss_weight", str(round(counts_loss_weight,2))]))
    file.write("\n")
    file.write("\t".join(["filters", str(args.filters)]))
    file.write("\n")
    file.write("\t".join(["n_dil_layers", str(args.n_dilation_layers)]))
    file.write("\n")
    file.write("\t".join(["bias_model_path", "{}bias_model_scaled.h5".format(args.output_prefix)]))
    file.write("\n")
    file.write("\t".join(["inputlen", str(args.inputlen)]))
    file.write("\n")
    file.write("\t".join(["outputlen", str(args.outputlen)]))
    file.write("\n")
    file.write("\t".join(["max_jitter", str(args.max_jitter)]))
    file.write("\n")
    file.write("\t".join(["chr_fold_path", str(args.chr_fold_path)]))
    file.write("\n")
    file.write("\t".join(["negative_sampling_ratio", str(args.negative_sampling_ratio)]))
    file.close()

if __name__=="__main__":
    # read the arguments
    parser = parse_data_args()
    args = parse_model_args(parser)

    main(args)

    
