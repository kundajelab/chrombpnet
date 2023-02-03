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
import json
import chrombpnet.training.utils.losses as losses
from chrombpnet.training.utils.data_utils import get_seq as get_seq
import chrombpnet.training.utils.one_hot as one_hot
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model


NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
PWM_SCHEMA = ["MOTIF_NAME", "MOTIF_PWM_FWD"]

def load_model_wrapper(args):
    # read .h5 model
    custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(args.model_h5)
    print("got the model")
    model.summary()
    return model


def fetch_footprinting_args():
    parser=argparse.ArgumentParser(description="get marginal footprinting for given model and given motifs")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Sequences and labels will be extracted centered at start (2nd col) + summit (10th col).")
    parser.add_argument("-fl", "--chr_fold_path", type=str, required=True, help="Path to file containing chromosome splits; we will only use the test chromosomes")
    parser.add_argument("-m", "--model_h5", type=str, required=True, help="Path to trained model, can be both bias or chrombpnet model")
    parser.add_argument("-bs", "--batch_size", type=int, default="64", help="input batch size for the model")
    parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
    parser.add_argument("-pwm_f", "--motifs_to_pwm", type=str, required=True, help="Path to a TSV file containing motifs in first column and motif string to use for footprinting in second column")    
    parser.add_argument("--ylim",default=None,type=tuple, required=False,help="lower and upper y-limits for plotting the motif footprint, in the form of a tuple i.e. \
    (0,0.8). If this is set to None, ylim will be autodetermined.")
    
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
    w_mot_seqs[:, midpoint-len(motif)//2:midpoint-len(motif)//2+len(motif)] = one_hot.dna_to_one_hot([motif])

    # midpoint of motif is the midpoint of sequence
    pred_output=model.predict(w_mot_seqs, batch_size=batch_size, verbose=True)
    footprint_for_motif_fwd = softmax(pred_output[0])*(np.exp(pred_output[1])-1)

    # reverse complement the sequence
    w_mot_seqs_revc = w_mot_seqs[:, ::-1, ::-1]
    pred_output_rev=model.predict(w_mot_seqs_revc, batch_size=batch_size, verbose=True)
    footprint_for_motif_rev = softmax(pred_output_rev[0])*(np.exp(pred_output_rev[1])-1)

    # add fwd sequence predictions and reverse sesquence predictions (not we flip the rev predictions)
    counts_for_motif = np.exp(pred_output_rev[1]) - 1 + np.exp(pred_output[1]) - 1
    footprint_for_motif_tot = footprint_for_motif_fwd+footprint_for_motif_rev[:,::-1]
    footprint_for_motif =  footprint_for_motif_tot / footprint_for_motif_tot.sum(axis=1)[:, np.newaxis]

    return footprint_for_motif.mean(0), counts_for_motif.mean(0)

def main(args):

	pwm_df = pd.read_csv(args.motifs_to_pwm, sep='\t',names=PWM_SCHEMA)
	print(pwm_df.head())
	genome_fasta = pyfaidx.Fasta(args.genome)

	model=load_model_wrapper(args)
	inputlen = model.input_shape[1] 
	outputlen = model.output_shape[0][1] 
	print("inferred model inputlen: ", inputlen)
	print("inferred model outputlen: ", outputlen)

	splits_dict = json.load(open(args.chr_fold_path))
	chroms_to_keep = set(splits_dict["test"])

	regions_df = pd.read_csv(args.regions, sep='\t', names=NARROWPEAK_SCHEMA)
	regions_subsample = regions_df[(regions_df["chr"].isin(chroms_to_keep))]
	regions_seqs = get_seq(regions_subsample, genome_fasta, inputlen)

	footprints_at_motifs = {}

	avg_response_at_tn5 = []

	motif = "control"
	motif_to_insert_fwd = ""
	print("inserting motif: ", motif)
	print(motif_to_insert_fwd)
	motif_footprint, motif_counts = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, model, inputlen, args.batch_size)
	footprints_at_motifs[motif]=[motif_footprint,motif_counts]

	plt.figure()
	plt.plot(range(200),motif_footprint[outputlen//2-100:outputlen//2+100])
	if args.ylim is not None: 
		plt.ylim(args.ylim)
	plt.xlabel("200bp arount motif insertion", fontsize=11)
	plt.ylabel("Probability", fontsize=11)
	plt.xticks(ticks=[0,100,200], labels=[-100,0,100])
	plt.tight_layout()
	plt.savefig(args.output_prefix+".{}.footprint.png".format(motif))


	#get motif names from column1 of the pwm_df
	for index, row in pwm_df.iterrows():
		motif=row["MOTIF_NAME"]
		motif_to_insert_fwd=row["MOTIF_PWM_FWD"]        
		print("inserting motif: ", motif)
		print(motif_to_insert_fwd)
		motif_footprint, motif_counts = get_footprint_for_motif(regions_seqs, motif_to_insert_fwd, model, inputlen, args.batch_size)
		footprints_at_motifs[motif]=[motif_footprint,motif_counts]

		# plot footprints of center 200bp
		if ("tn5" in motif.lower()) or ("dnase" in motif.lower()):
				avg_response_at_tn5.append(np.round(np.max(motif_footprint),3))
		plt.figure()
		plt.plot(range(200),motif_footprint[outputlen//2-100:outputlen//2+100])
		if args.ylim is not None: 
			plt.ylim(args.ylim)
		plt.xlabel("200bp arount motif insertion", fontsize=11)
		plt.ylabel("Probability", fontsize=11)
		plt.xticks(ticks=[0,100,200], labels=[-100,0,100])
		plt.tight_layout()
		plt.savefig(args.output_prefix+".{}.footprint.png".format(motif))

	if len(avg_response_at_tn5) > 0:
		if np.all(np.array(avg_response_at_tn5) < 0.003):
			ofile = open("{}_max_bias_resonse.txt".format(args.output_prefix), "w")
			ofile.write("corrected_"+str(round(np.mean(avg_response_at_tn5),3))+"_"+"/".join(list(map(str,avg_response_at_tn5))))
			ofile.close()
		else:
			ofile = open("{}_max_bias_resonse.txt".format(args.output_prefix), "w")
			ofile.write("uncorrected_"+str(round(np.mean(avg_response_at_tn5),3))+"_"+"/".join(list(map(str,avg_response_at_tn5))))
			ofile.close()

	print("Saving marginal footprints")
	dd.io.save("{}_footprints.h5".format(args.output_prefix),
		footprints_at_motifs,
		compression='blosc')


if __name__ == '__main__':
    args=fetch_footprinting_args()
    main(args)
