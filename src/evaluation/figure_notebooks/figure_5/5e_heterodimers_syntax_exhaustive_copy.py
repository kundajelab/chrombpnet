#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"]="3"
#get_ipython().run_line_magic('matplotlib', 'inline')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42 


# In[2]:


#regions = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/negatives_data/negatives_with_summit.bed"
#regions = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/IMR90/negatives_data/negatives_with_summit.bed"
#regions = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/negatives_data/negatives_with_summit.bed"
regions = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/IMR90/negatives_data/negatives_with_summit.bed"

genome = "/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa"
#model_h5 = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/K562/nautilus_runs_may18/K562_05.13.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5"
#model_h5 = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/K562/nautilus_runs/K562_02.17.2022_bias_128_4_1234_0.5_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5"
#model_h5="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/IMR90/nautilus_runs_apr12/IMR90_04.09.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5"
model_h5="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/IMR90/nautilus_runs_apr12/IMR90_04.09.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5"
#model_h5 ="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/GM12878/nautilus_runs/GM12878_03.06.2022_bias_128_4_1234_0.8_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5"
#model_h5 = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/nautilus_runs/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/chrombpnet_model/chrombpnet_wo_bias.h5"


# In[3]:


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
    pred_output=model.predict(w_mot_seqs, batch_size=batch_size, verbose=False)
    footprint_for_motif_fwd = softmax(pred_output[0])*(np.exp(pred_output[1])-1)

    # reverse complement the sequence
    w_mot_seqs_revc = w_mot_seqs[:, ::-1, ::-1]
    pred_output_rev=model.predict(w_mot_seqs_revc, batch_size=batch_size, verbose=False)
    footprint_for_motif_rev = softmax(pred_output_rev[0])*(np.exp(pred_output_rev[1])-1)

    # add fwd sequence predictions and reverse sesquence predictions (not we flip the rev predictions)
    counts_for_motif = (np.exp(pred_output_rev[1]) - 1 + np.exp(pred_output[1]) - 1)/2
    footprint_for_motif_tot = (footprint_for_motif_fwd+footprint_for_motif_rev[:,::-1])/2
    footprint_for_motif =  footprint_for_motif_tot / footprint_for_motif_tot.sum(axis=1)[:, np.newaxis]

    return footprint_for_motif_tot, counts_for_motif 

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
    
    #print(motif,start,start+len(motif))
    
    # midpoint of motif is the midpoint of sequence
    pred_output=model.predict(w_mot_seqs, batch_size=batch_size, verbose=False)
    footprint_for_motif_fwd = softmax(pred_output[0])*(np.exp(pred_output[1])-1)

    # reverse complement the sequence
    w_mot_seqs_revc = w_mot_seqs[:, ::-1, ::-1]
    pred_output_rev=model.predict(w_mot_seqs_revc, batch_size=batch_size, verbose=False)
    footprint_for_motif_rev = softmax(pred_output_rev[0])*(np.exp(pred_output_rev[1])-1)

    # add fwd sequence predictions and reverse sesquence predictions (not we flip the rev predictions)
    counts_for_motif = (np.exp(pred_output_rev[1]) - 1 + np.exp(pred_output[1]) - 1)/2
    footprint_for_motif_tot = (footprint_for_motif_fwd+footprint_for_motif_rev[:,::-1])/2
    footprint_for_motif =  footprint_for_motif_tot / footprint_for_motif_tot.sum(axis=1)[:, np.newaxis]

    return footprint_for_motif_tot, counts_for_motif


def new_orient_get_footprint_for_two_motifs(seqs, motifs, model, inputlen, batch_size, spacing):
    '''
    Returns footprints for a given motif. Motif is inserted in both the actual sequence and reverse complemented version.
    seqs input is already assumed to be one-hot encoded. motif is in sequence format.
    '''
    midpoint=inputlen//2

    spacing_per_motif = spacing // 2
    
    w_mot_seqs = seqs.copy()
    
    motif = motifs[0]
    start = midpoint-(len(motif)//2)
    w_mot_seqs[:, start:start+len(motif)] = dinuc_shuffle_main.dna_to_one_hot([motif])[:, ::-1, ::-1]
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
    
    #print(motif,start,start+len(motif))
    
    # midpoint of motif is the midpoint of sequence
    pred_output=model.predict(w_mot_seqs, batch_size=batch_size, verbose=False)
    footprint_for_motif_fwd = softmax(pred_output[0])*(np.exp(pred_output[1])-1)

    # reverse complement the sequence
    w_mot_seqs_revc = w_mot_seqs[:, ::-1, ::-1]
    pred_output_rev=model.predict(w_mot_seqs_revc, batch_size=batch_size, verbose=False)
    footprint_for_motif_rev = softmax(pred_output_rev[0])*(np.exp(pred_output_rev[1])-1)

    # add fwd sequence predictions and reverse sesquence predictions (not we flip the rev predictions)
    counts_for_motif = (np.exp(pred_output_rev[1]) - 1 + np.exp(pred_output[1]) - 1)/2
    footprint_for_motif_tot = (footprint_for_motif_fwd+footprint_for_motif_rev[:,::-1])/2
    footprint_for_motif =  footprint_for_motif_tot / footprint_for_motif_tot.sum(axis=1)[:, np.newaxis]

    return footprint_for_motif_tot, counts_for_motif


# In[4]:


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


# In[5]:


model=load_model(model_h5)


# In[6]:


NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
inputlen = 2114
regions_df = pd.read_csv(regions, sep='\t', names=NARROWPEAK_SCHEMA)
chroms_to_keep = ["chr1"]
regions_subsample = regions_df[(regions_df["chr"].isin(chroms_to_keep))].sample(100, random_state=0)
genome_fasta = pyfaidx.Fasta(genome)
regions_seqs = get_seq(regions_subsample, genome_fasta, inputlen)


# In[7]:


motif =  ""
batch_size=128
full_footprint_control = get_footprint_for_motif(regions_seqs, motif, model, inputlen, batch_size)


# In[8]:


batch_size=512

def simulate_motif_set(motifs, model, regions_seqs):

    
    #print(motifs[0])
    #print(motifs[1])
    #motif1_only = get_footprint_for_motif(regions_seqs, motifs[0], model, inputlen, batch_size)[1]
    #motif2_only = get_footprint_for_motif(regions_seqs, motifs[1], model, inputlen, batch_size)[1]
    #total_sum = motif1_only+motif2_only
    
    #plt.plot(motif1_only)
    #plt.plot(motif2_only)

    # -/- orientation
    motif1_only = []
    motif2_only = []
    data_in_spacings = []
    for spacing in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,40,50,60,70,80,90,100,110,120,130,140,150]:
    #for spacing in [4,6]:

        puu_runx_footprint = get_footprint_for_two_motifs(regions_seqs, [motifs[0], ""], model, inputlen, batch_size, spacing=-1*spacing)
        motif1_only.append(puu_runx_footprint)

        puu_runx_footprint = get_footprint_for_two_motifs(regions_seqs, ["", motifs[1]], model, inputlen, batch_size, spacing=-1*spacing)
        motif2_only.append(puu_runx_footprint)

        puu_runx_footprint = get_footprint_for_two_motifs(regions_seqs, motifs, model, inputlen, batch_size, spacing=-1*spacing)
        data_in_spacings.append(puu_runx_footprint)
        
    # +/+ orientation
    motif1_only_rev = []
    motif2_only_rev = []
    data_in_spacings_rev = []
    for spacing in [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,40,50,60,70,80,90,100,110,120,130,140,150]:
    #for spacing in [4,6]:
        puu_runx_footprint = get_footprint_for_two_motifs(regions_seqs, [motifs[0], ""], model, inputlen, batch_size, spacing=spacing)
        motif1_only_rev.append(puu_runx_footprint)

        puu_runx_footprint = get_footprint_for_two_motifs(regions_seqs, ["", motifs[1]], model, inputlen, batch_size, spacing=spacing)
        motif2_only_rev.append(puu_runx_footprint)
    
        puu_runx_footprint = get_footprint_for_two_motifs(regions_seqs, motifs, model, inputlen, batch_size, spacing=spacing)
        data_in_spacings_rev.append(puu_runx_footprint)
        
        
    #all_counts = [x[1] for x in data_in_spacings]
    #all_counts_rev = [x[1] for x in data_in_spacings_rev]

    distance = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,40,50,60,70,80,90,100,110,120,130,140,150]+[-1,-2,-3,-4,-5,-6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-20,-40,-50,-60,-70,-80,-90,-100,-110,-120,-130,-140,-150]

    #total_counts = all_counts_rev+all_counts
    #distance = np.array(distance)
    #total_counts = np.array(total_counts)
    
    #coop_score = np.max(total_counts)
    #coop_distance = distance[np.argmax(total_counts)]
   
    #print(motif1_only, motif2_only, total_sum, coop_score, np.max(total_counts)/np.mean(total_counts), coop_distance, total_counts,distance)

    return [motif1_only, motif2_only, data_in_spacings, motif1_only_rev, motif2_only_rev, data_in_spacings_rev, distance]

    #print(motif1_only, motif2_only, total_sum, coop_score, np.max(total_counts)/np.mean(total_counts), coop_distance)

    


# In[9]:


import pandas as pd


# In[10]:


annotations = pd.read_csv("addons/imr90/imr90.counts.tomtom.motifs_string.tsv",header=0, sep="\t")


# In[11]:


annotations.head()


# In[13]:


def define_cooperative_net_pos_neg(motif1_only, motif2_only, cooop):
    total_sum_th = motif1_only + motif2_only
    total_diff_th = np.abs(motif1_only - motif2_only)
    total_diff_th1 = np.max([motif1_only,motif2_only])
    total_sum_th1 = np.max([motif1_only,motif2_only])
    
    pos_coop_score = (np.max(cooop) - total_sum_th)/total_sum_th
    pos_coop_score1 = (np.max(cooop) - total_sum_th1)/total_sum_th1

    neg_coop_score = (np.min(cooop) - total_diff_th)/total_diff_th
    neg_coop_score1 = (np.min(cooop) - total_diff_th1)/total_diff_th1
        
    return pos_coop_score, pos_coop_score1, neg_coop_score, neg_coop_score1


# In[ ]:


#motifs = ["TGGAC","TGACTCAT"]
#motifs = ["AGGAATGT","TTGACTCA"]
values_fresh = {}

#for i in range(annotations.shape[0]-1):
#    for j in range(i+1,annotations.shape[0]):
for i in range(annotations.shape[0]):
    for j in range(i,annotations.shape[0]):
        label1 = annotations.loc[i,"Pattern"] + "-" + annotations.loc[j,"Pattern"]
        label2 = annotations.loc[i,"Label"] + "-" + annotations.loc[j,"Label"]
        motifs = [annotations.loc[i,"string"],annotations.loc[j,"string"]]

        #print(label1)
        #print(label2)
        #print(motifs)
        #list_motifs = [["CAGGTG","TGACTCAT"],["AGGAATGT","TTGACTCA"]]
        #labels = ["ZEB-AP1", "TEAD-AP1"]
        values_fresh[label1] = [simulate_motif_set(motifs, model, regions_seqs),label2]
        #break
   #break


# In[ ]:


motif =  ""
batch_size=128
full_footprint_control = get_footprint_for_motif(regions_seqs, motif, model, inputlen, batch_size)
values_fresh["control"] = full_footprint_control


# In[ ]:

import pickle
with open('addons/imr90/coop_matrix_all_new_spaced_3.pkl', 'wb') as handle:
    pickle.dump(values_fresh, handle)

