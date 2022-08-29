import pyBigWig
import pandas as pd
import numpy as np
import deepdish as dd
import os
import pyfaidx
import random
import pickle as pkl
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import tensorflow as tf
import argparse
import context
import json
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import tensorflow as tf
import deepdish

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
PWM_SCHEMA = ["MOTIF_NAME", "MOTIF_PWM_FWD", "1", "2", "3", "4", "5", "6", "7", "8"]


def load_model_wrapper(args):
    # read .h5 model
    custom_objects={ "tf": tf}    
    get_custom_objects().update(custom_objects)    
    corrected_model=load_model(args.corrected_model_h5, compile=False)
    uncorrected_model=load_model(args.uncorrected_model_h5, compile=False)
    print("got the model")
    corrected_model.summary()
    return corrected_model,uncorrected_model


def get_footprint_characteristics(motif_width, corrected_footprint_1, nomotif_corrected_footprint_1):


    fc_change = corrected_footprint_1 -  nomotif_corrected_footprint_1
    smoothed_fc = np.convolve(fc_change, np.ones(5)/5, mode='same')
    grad = np.gradient(smoothed_fc)

    rs = 500-motif_width//2-50
    ls = 500+motif_width//2+50

    midpoint=motif_width//2+50
    ldx = midpoint+np.argmax(grad[500:ls])
    rdx = np.argmin(grad[rs:500])
    minv = np.min(smoothed_fc[rs:ls])
    maxv = np.max(smoothed_fc[rs:ls])

    #offset_right = np.min((grad[rs:ls][0:rdx][::-1]>0).nonzero())
    #offset_left = np.min((grad[rs:ls][ldx:]<0).nonzero())

    #right_max = rdx - offset_right
    #left_max = ldx + offset_left

    #maxv = np.max([smoothed_fc[rs:ls][right_max], smoothed_fc[rs:ls][left_max]])

    #scores_path = "{}_{}_footprint_score.txt".format(output_prefix, motif_name)

    #f =  open(scores_path, "w")
    text = ",".join([str(np.round(maxv-minv,2))])

    #f.write(text)
    #f.close()

    return text

def fetch_footprinting_args():
    parser=argparse.ArgumentParser(description="get marginal footprinting for given model and given motifs")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Path to file containing chromosome splits; we will only use the test chromosomes")
    parser.add_argument("-cm", "--corrected_model_h5", type=str, required=True, help="Path to trained model, can be both bias or chrombpnet model")
    parser.add_argument("-um", "--uncorrected_model_h5", type=str, required=True, help="Path to trained model, can be both bias or chrombpnet model")
    parser.add_argument("-bs", "--batch_size", type=int, default="64", help="input batch size for the model")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-tt", "--title", type=str, required=True, help="title for plot")
    parser.add_argument("-pwm_f", "--motifs_to_pwm", type=str, required=True, 
                        help="Path to a TSV file containing motifs in first column and motif string to use for footprinting in second column")    
    
    args = parser.parse_args()
    return args

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def get_footprint_for_motif(seqs, motif, model, inputlen, batch_size):
    '''
    Returns footprints for a given motif. Motif is inserted in both the actual sequence and reverse complemented version.
    seqs input is already assumed to be one-hot encoded. motif is in sequence format.
    '''
    midpoint=inputlen//2

    w_mot_seqs = seqs.copy()
    w_mot_seqs[:, midpoint-len(motif)//2:midpoint-len(motif)//2+len(motif)] = context.one_hot.dna_to_one_hot([motif])

    # midpoint of motif is the midpoint of sequence
    pred_output=model.predict(w_mot_seqs, batch_size=batch_size, verbose=True)
    footprint_for_motif_fwd = softmax(pred_output[0])*(np.exp(pred_output[1])-1)

    # reverse complement the sequence
    w_mot_seqs_revc = w_mot_seqs[:, ::-1, ::-1]
    pred_output_rev=model.predict(w_mot_seqs_revc, batch_size=batch_size, verbose=True)
    footprint_for_motif_rev = softmax(pred_output_rev[0])*(np.exp(pred_output_rev[1])-1)

    # add fwd sequence predictions and reverse sesquence predictions (not we flip the rev predictions)
    counts_for_motif = (np.exp(pred_output_rev[1]) - 1 + np.exp(pred_output[1]) - 1)/2
    footprint_for_motif_tot = (footprint_for_motif_fwd+footprint_for_motif_rev[:,::-1])/2
    footprint_for_motif =  footprint_for_motif_tot / footprint_for_motif_tot.sum(axis=1)[:, np.newaxis]

    return footprint_for_motif_tot.mean(0), counts_for_motif.mean(0), footprint_for_motif_tot

def main():

    args=fetch_footprinting_args()

    pwm_df = pd.read_csv(args.motifs_to_pwm, sep='\t',names=PWM_SCHEMA)
    print(pwm_df)
    genome_fasta = pyfaidx.Fasta(args.genome)

    model,uncorrected_model=load_model_wrapper(args)
    inputlen = model.input_shape[1] 
    outputlen = model.output_shape[0][1] 
    print("inferred model inputlen: ", inputlen)
    print("inferred model outputlen: ", outputlen)

    splits_dict = json.load(open(args.chr_fold_path))
    chroms_to_keep = set(splits_dict["test"])

    regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
    chroms_to_keep = ["chr1"]
    regions_subsample = regions_df[(regions_df["chr"].isin(chroms_to_keep))].sample(1000)
    regions_seqs = context.get_seq(regions_subsample, genome_fasta, inputlen)

    footprints_at_motifs = {}
    footprints_at_motifs_uncorrected = {}

    # get control sequence
    motif="control"
    motif_to_insert_fwd = ""

    motif_footprint, motif_counts, motif_vals = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, model, inputlen, args.batch_size)
    uncorrecetd_motif_footprint, uncorrected_motif_counts, uncorrected_motif_vals = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, uncorrected_model, inputlen, args.batch_size)

    footprints_at_motifs[motif]=[motif_footprint,motif_counts, motif_vals]
    #footprints_at_motifs_uncorrected[motif]=[uncorrecetd_motif_footprint,uncorrected_motif_counts, uncorrected_motif_vals]
    values = []
    for x in footprints_at_motifs[motif][2]:
        val = get_footprint_characteristics(0, x, footprints_at_motifs[motif][0])
        values.append(float(val.strip()))

    print(np.quantile(values, 0.95))
        


if __name__ == '__main__':
    main()

