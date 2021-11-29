
import pyBigWig
import pandas as pd
import numpy as np
import deepdish as dd
import os
import pyfaidx
import random
import pickle as pkl
import matplotlib.pyplot as plt
import tensorflow as tf
import argparse
import context


NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
PWM_SCHEMA = ["MOTIF_NAME", "MOTIF_PWM_FWD"]


def fetch_footprinting_args():
    parser=argparse.ArgumentParser(description="get marginal footprinting for given model and given motifs")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-chr", "--test_chr", type=str, required=True, help="Only use the following chromsomes from the regions argument")
    parser.add_argument("-m", "--model_h5", type=str, required=True, help="Path to trained model, can be both bias or chrombpnet model")
    parser.add_argument("-bs", "--batch_size", type=int, default="64", help="input batch size for the model")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-pwm_f", "--motifs_to_pwm", type=str, required=True, 
                        help="Path to a TSV file containing motifs in first column and motif string to use for footprinting in second column")    
    parser.add_argument("-mo", "--motifs", nargs="+", type=lambda s: [str(item) for item in s.split(',')], default=["Tn5"],
                        help="Input motifs to do marginal footprinting, the motif names input here should be present in the collumn one in motifs_to_pwm file argument")
    
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

    #midpoint of motif is the midpoint of sequence
    output=model.predict(w_mot_seqs, batch_size=batch_size, verbose=True)[0]
    footprint_for_motif_fwd = softmax(output)

    # reverse complement the sequence
    w_mot_seqs_revc = w_mot_seqs[:, ::-1, ::-1]
    footprint_for_motif_rev = softmax(model.predict(w_mot_seqs_revc, batch_size=batch_size, verbose=True)[0])  

    # add fwd sequence predictions and reverse sesquence predictions (not we flip the rev predictions)
    footprint_for_motif = footprint_for_motif_fwd+footprint_for_motif_rev[:,::-1,:]

    return footprint_for_motif.mean(0)

def main():

    args=fetch_footprinting_args()

    pwm_df = pd.read_csv(args.motifs_to_pwm, sep='\t',names=PWM_SCHEMA)
    print(pwm_df)
    genome_fasta = pyfaidx.Fasta(args.genome)

    model=context.load_model_wrapper(args)
    inputlen = model.input_shape[1] 
    outputlen = model.output_shape[0][1] 
    print("inferred model inputlen: ", inputlen)
    print("inferred model outputlen: ", outputlen)

    regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
    regions_subsample = regions_df[(regions_df["chr"]==args.test_chr)]
    regions_seqs = context.get_seq(regions_subsample, genome_fasta, inputlen)

    footprints_at_motifs = {}

    for motif in args.motifs[0]:
        print("inserting motif: ", motif)
        motif_to_insert_fwd = pwm_df[pwm_df["MOTIF_NAME"]==motif]["MOTIF_PWM_FWD"].values[0]
        print(motif_to_insert_fwd)
        motif_footprint = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, model, inputlen, args.batch_size)
        footprints_at_motifs[motif]=motif_footprint

        # plot footprints of center 200bp
        plt.figure()
        plt.plot(range(200),motif_footprint[outputlen//2-100:outputlen//2+100])
        plt.savefig(args.output_prefix+".{}.footprint.png".format(motif))

    print("Saving marginal footprints")
    dd.io.save("{}.footprints.h5".format(args.output_prefix),
        footprints_at_motifs,
        compression='blosc')



if __name__ == '__main__':
    main()

