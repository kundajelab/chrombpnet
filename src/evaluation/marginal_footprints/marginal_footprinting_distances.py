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
import one_hot as dinuc_shuffle_main
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import tensorflow as tf
import os
os.environ["CUDA_VISIBLE_DEVICES"]="1"

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
    w_mot_seqs[:, midpoint-len(motif)//2:midpoint-len(motif)//2+len(motif)] =dinuc_shuffle_main.dna_to_one_hot([motif])

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
    k_l=1000
    # if counts_for_motif.shape[0] > k_l:
    #     print(counts_for_motif.shape)
    #     smallest_100 = np.argpartition(counts_for_motif[:,0], k_l)[:k_l]
    #     return footprint_for_motif_tot[smallest_100].mean(0), counts_for_motif[smallest_100].mean(0), smallest_100
    # else:
    return footprint_for_motif.mean(0), counts_for_motif.mean(0), None

def get_footprint_for_two_motifs(seqs, motifs, model, inputlen, batch_size, spacing):
    '''
    Returns footprints for a given motif. Motif is inserted in both the actual sequence and reverse complemented version.
    seqs input is already assumed to be one-hot encoded. motif is in sequence format.
    '''
    midpoint=inputlen//2

    spacing_per_motif = spacing // 2
    
    w_mot_seqs = seqs.copy()
    
    motif = motifs[0]
    start = midpoint-(len(motif)//2)
    w_mot_seqs[:, start:start+len(motif)] = dinuc_shuffle_main.dna_to_one_hot([motif])
    #print(motif,start,start+len(motif))
    if spacing > 0:
        spacing_per_motif = spacing 
        motif = motifs[1]
        start = start+len(motifs[0])+spacing_per_motif 
        w_mot_seqs[:, start:start+len(motif)] = dinuc_shuffle_main.dna_to_one_hot([motif])
    else:
        spacing_per_motif = spacing 
        motif = motifs[1]
        start = start + spacing_per_motif - len(motif)
        w_mot_seqs[:, start:start+len(motif)] = dinuc_shuffle_main.dna_to_one_hot([motif])
    

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

    return footprint_for_motif.mean(0), counts_for_motif.mean(0)

def get_seq(peaks_df, genome, width, shuffle=False):
    """
    fetches sequence from a given genome.
    """
    vals = []

    for i, r in peaks_df.iterrows():
        sequence = str(genome[r['chr']][(r['start']+r['summit'] - width//2):(r['start'] + r['summit'] + width//2)])
        if len(sequence) == width:
                vals.append(sequence)

    return dinuc_shuffle_main.dna_to_one_hot(vals)

genome = "/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa"
genome_fasta = pyfaidx.Fasta(genome)
NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
inputlen = 2114

modelcsv = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/v1/model_dir_atac.csv", sep=",", header=None)
motif_table = pd.read_csv("composites_summary_table_2.tsv", sep='\t', header=0)

store_data_object = {}
for i,r in motif_table.iterrows():
	celltype = r['model'].split("_")[0]
	#print(celltype)
	model_dir=modelcsv[(modelcsv[1]==celltype) & (modelcsv[0]=="fold_0")][2].values[0]
	negstr="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/celltype/negatives_data/negatives_with_summit.bed"
	negdir=negstr.replace("celltype", celltype)
	regions_df = pd.read_csv(negdir, sep='\t', names=NARROWPEAK_SCHEMA)
	mdl_path1 = os.path.join(model_dir,"chrombpnet_model/chrombpnet_wo_bias.h5")
	model=load_model(mdl_path1, compile=False)

	data_in_spacings = {}
	counts_in_spacings = {}

	regions_subsample = regions_df[regions_df['chr']=="chr1"].sample(100, random_state=0)
	regions_seqs = get_seq(regions_subsample, genome_fasta, inputlen)

	motif =  ""
	batch_size=128
	control_footprint = get_footprint_for_motif(regions_seqs, motif, model, inputlen, batch_size)

	motif=r['motif1_strings_corrected']
	batch_size=128
	motif1_footprint = get_footprint_for_motif(regions_seqs, motif, model, inputlen, batch_size)

	motif=r['motif2_strings_corrected']
	batch_size=128
	motif2_footprint = get_footprint_for_motif(regions_seqs, motif, model, inputlen, batch_size)

	motifs = [r['motif1_strings_corrected'], r['motif2_strings_corrected']]

	for spacing in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,40,50,75,100,125, 150]:
		full_footprint_data = get_footprint_for_two_motifs(regions_seqs, motifs, model, inputlen, batch_size, spacing=-1*spacing)
		data_in_spacings[str(-1*spacing)] = full_footprint_data[0]
		counts_in_spacings[str(-1*spacing)] = full_footprint_data[1]


	for spacing in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,40,50,75,100,125, 150]:
		full_footprint_data = get_footprint_for_two_motifs(regions_seqs, motifs, model, inputlen, batch_size, spacing=spacing)
		data_in_spacings[str(spacing)] = full_footprint_data[0]
		counts_in_spacings[str(spacing)] = full_footprint_data[1]

	#print("num distances",len(data_in_spacings.keys()))
	#print(data_in_spacings.keys())

	store_data_object[r['name_x']] = {'control':{'counts': control_footprint[1], 'profile':control_footprint[0]}, 
							'motif1':{'counts': motif1_footprint[1], 'profile': motif1_footprint[0] }, 
							'motif2': {'counts': motif2_footprint[1], 'profile': motif2_footprint[0]},
						'combine_footprint':data_in_spacings, 'combine_counts':counts_in_spacings	}

	#print(counts_in_spacings['5'], r['name_x'])
	#break 
	#if i==2:
	#	break

import json
import h5py

# Save function
def save_dict_to_h5(h5group, dic):
    for key, item in dic.items():
        if isinstance(item, dict):
            subgroup = h5group.create_group(key)
            save_dict_to_h5(subgroup, item)
        elif isinstance(item, list):
            subgroup = h5group.create_group(key)
            for idx, sub_item in enumerate(item):
                subgroup.create_dataset(str(idx), data=sub_item)
        else:
            h5group.create_dataset(key, data=item)


# Function to save the whole dictionary to HDF5 file
def save_to_h5(filename, data):
    with h5py.File(filename, 'w') as h5file:
        save_dict_to_h5(h5file, data)


save_to_h5('composite_footprint_data.h5', store_data_object)









