import pyfaidx
import math
import pandas as pd
import numpy as np
import pyBigWig
from modisco.visualization import viz_sequence
import matplotlib.pyplot as plt
import argparse
import one_hot
import os

def parse_args():
    parser=argparse.ArgumentParser(description="build pwm matrix from bigwig")
    parser.add_argument("-i","--bigwig", help="generated bigiwig file")
    parser.add_argument("-b","--bed", help="build pwm from the regions in this bed, regions are centered at the midpoint of start,end")
    parser.add_argument("-g", "--ref_fasta", help="reference genome fasta")
    parser.add_argument("-o", "--out_dir", help="output dir for storing pwm")
    parser.add_argument("-c","--chr",type=str, default="chr20")
    return parser.parse_args()

args = parse_args()

hg38 = pyfaidx.Fasta(args.ref_fasta)
bw = pyBigWig.open(args.bigwig) 
gc_neg_peaks =  pd.read_csv(args.bed,
                           sep='\t',
                           names=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"])
gc_neg_peaks.head()

peaks = gc_neg_peaks[(gc_neg_peaks["chr"]==args.chr)]

def get_seq(genome, peaks_df, width=2000):
    vals = []
    for i, r in peaks_df.iterrows():
        seq = str(genome[r['chr']][(r['start'] + r['summit']) - width//2:(r['start'] + r['summit']) + width//2])
        if len(seq) == 2000:
            vals.append(seq)
        else:
            #print(len(seq))
            #print((r['start'] + r['summit'])- width//2, (r['start'] + r['summit']) + width//2)
            pass
        
    return one_hot.dna_to_one_hot(vals)

def get_cts(bw, peaks_df, width=2000):
    vals = []
    for i, r in peaks_df.iterrows():
        try:
            vals.append(np.nan_to_num(bw.values(r['chr'], 
                                                (r['start'] + r['summit'])- width//2,
                                                (r['start'] + r['summit']) + width//2)))
        except:
            #print(i)
            #print((r['start'] + r['summit'])- width//2, (r['start'] + r['summit']) + width//2)
            pass
        
    return np.array(vals)

seqs = get_seq(hg38, peaks)
cnts = get_cts(bw, peaks)

def get_pwm_bg(seqs, cnts):
    new_seqs = []
    for i in range(cnts.shape[0]):
        for j in range(13,cnts.shape[1]-13):
            if cnts[i,j] > 0:
                new_seqs.append(seqs[i,j-12:j+12,:]* cnts[i,j])
    #print(np.array(new_seqs).shape)
    motif = np.sum(new_seqs, axis=0)
    motif = motif/np.sum(motif, axis=-1, keepdims=True)
    bg = np.sum(np.sum(new_seqs, axis=0), axis=0)
    bg = bg/sum(bg)
    return motif, bg

motif, bg = get_pwm_bg(seqs, cnts)


figsize=(20,2)
fig = plt.figure(figsize=figsize)
ax = fig.add_subplot(111) 
viz_sequence.plot_weights_given_ax(ax=ax, array=viz_sequence.ic_scale(motif, background=bg))
plt.savefig(os.path.join(args.out_dir, "bias_pwm.png"))
