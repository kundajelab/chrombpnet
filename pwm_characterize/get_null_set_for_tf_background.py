#This script does the following:
#rank peaks by observed count signal
#get the peaks with low count signal (i.e. bottom 5%)
#dinuc shuffle those bottom 5%
#score the dinuc shuffled seqs with model
#compute entropy
#rank by entropy, take 100 with highest entropy

import argparse
import pandas as pd
import numpy as np 
import pyBigWig
import pysam
import pickle

#import the libraries necessary for getting keras model predictions
from scipy.special import softmax,expit
from scipy.stats import entropy 
from kerasAC.interpret.deepshap import * 
from kerasAC.interpret.profile_shap import * 
from kerasAC.vis import * 
from kerasAC.helpers.transform_bpnet_io import * 
from kerasAC.util import * 
#load the model! 
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import * 
from kerasAC.custom_losses import *

def load_model_wrapper(model_path):
    custom_objects={"recall":recall,
                    "sensitivity":recall,
                    "specificity":specificity,
                    "fpr":fpr,
                    "fnr":fnr,
                    "precision":precision,
                    "f1":f1,
                    "ambig_binary_crossentropy":ambig_binary_crossentropy,
                    "ambig_mean_absolute_error":ambig_mean_absolute_error,
                    "ambig_mean_squared_error":ambig_mean_squared_error,
                    "MultichannelMultinomialNLL":MultichannelMultinomialNLL}
    get_custom_objects().update(custom_objects)
    model=load_model(model_path)
    return model

def get_preds(model,seq_onehot):
    preds=model.predict(seq_onehot)
    prof=preds[0] #logits
    probs=np.squeeze(softmax(prof,axis=1)) #probabilities 
    count=np.squeeze(preds[1])   #count (1 val) 
    count_track=probs*np.expand_dims(count,axis=1)  #count track
    return probs,count_track

def dinuc_shuffle(seq):
    #get list of dinucleotides
    nucs=[]
    for i in range(0,len(seq),2):
        nucs.append(seq[i:i+2])
    #generate a random permutation
    #set the seed so this shuffling is reproducible for a given sequence
    random.seed(1234) 
    random.shuffle(nucs)
    return ''.join(nucs) 



def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('--peak_file')
    parser.add_argument('--bigwig_count_file')
    parser.add_argument("--ref_fasta",default="/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta")
    parser.add_argument('--bpnet_out_flank',type=int,default=500)
    parser.add_argument('--bpnet_input_flank',type=int,default=1057)
    parser.add_argument("--model_path")
    parser.add_argument('--n_to_sample',type=int,default=100)
    parser.add_argument('--out_pickle')
    return parser.parse_args()

def main():
    args=parse_args()
    #load the keras model
    model=load_model_wrapper(args.model_path)
    print("loaded model") 
    #open reference fasta file 
    ref=pysam.FastaFile(args.ref_fasta)
    #open bigwig signal track for reading
    signal=pyBigWig.open(args.bigwig_count_file,'r')
    #read in the peak file
    peaks=pd.read_csv(args.peak_file,header=None,sep='\t')
    print("loaded peak file")
    peak_to_count_signal={}
    #get the output region counts for each peak
    peaks['summit_pos']=peaks[1]+peaks[9]
    peaks['start_pos']=peaks['summit_pos']-args.bpnet_out_flank
    peaks['end_pos']=peaks['summit_pos']+args.bpnet_out_flank
    print("calculating the total count signal in the peaks") 
    for index,row in peaks.iterrows():
        cur_signal=np.sum(signal.values(row[0],row['start_pos'],row['end_pos']))
        peak_to_count_signal[(row[0],row[1],row[2],row[9])]=cur_signal
    print("got the count signal for all the peak tracks")
    #sort the peaks by count signal from lowest to highest, and select the bottom 5%
    sorted_peak_to_count=sorted(peak_to_count_signal.items(), key=lambda item: item[1])
    num_peaks=len(sorted_peak_to_count)
    #get the 5% of peaks with the lowest signal
    lowest_signal_n=round(0.05*num_peaks)
    peaks_with_low_signal=sorted_peak_to_count[0:lowest_signal_n]
    #dinucleotide shuffle these peaks, keep the same order as the "peaks_with_low_signal" list so that we can map back directly  
    shuffled_seqs={}
    peak_index=0 
    for peak in peaks_with_low_signal:
        if peak_index%100==0:
            print(str(peak_index)+'/'+str(peaks_with_low_signal))
        peak_info=peak[0]
        count_signal=peak[1]
        chrom=peak_info[0]
        summit=peak_info[1]+peak_info[3] 
        input_seq_start=summit-args.bpnet_input_flank
        input_seq_end=summit+args.bpnet_input_flank 
        ref_seq=ref.fetch(chrom,input_seq_start,input_seq_end)
        shuffled_seq=dinuc_shuffle(ref_seq)
        shuffled_seq_onehot=one_hot_encode([seq])
        probability_preds,count_preds=get_preds(model,shuffled_seq_onehot)
        #get the entropy from the probability track
        signal_entropy=entropy(probability_preds) 
        #store the count track so we can subtract it out as baseline when we insert the motifs into the dinuc shuffled sequence 
        shuffled_seq[(peak_info[0],peak_info[1],peak_info[2],peak_info[3],count_signal,ref_seq,shuffled_seq,count_preds,probability_preds)]=signal_entropy
        
    #sort the peaks by entropy, store the n peaks with the highest entropy
    shuffled_seq_entropy_sorted=sorted(shuffled_seq.items(), key=lambda item: item[1],decreasing=True)
    selected_background=shuffled_seq_entropy_sorted[0:args.n_to_sample]
    #store the background as a dictionary
    print("prepping the outputs") 
    background_dict={}
    for entry in selected_background:
        seq_info=entry[0]
        entropy=entry[1]
        chrom=seq_info[0]
        start_pos=seq_info[1]
        end_pos=seq_info[2]
        summit=seq_info[3]
        background_dict[(chrom,start_pos,end_pos,summit)]={}
        background_dict[(chrom,start_pos,end_pos,summit)]['count_signal']=seq_info[4]
        background_dict[(chrom,start_pos,end_pos,summit)]['ref_seq']=seq_info[5]
        background_dict[(chrom,start_pos,end_pos,summit)]['shuffled_seq']=seq_info[6]
        background_dict[(chrom,start_pos,end_pos,summit)]['count_preds']=seq_info[7]
        background_dict[(chrom,start_pos,end_pos,summit)]['probability_preds']=seq_info[8]
        background_dict[(chrom,start_pos,end_pos,summit)]['signal_entropy']=entropy
    
    #pickle the background dict!
    print("generating the pickle") 
    with open(args.output_pickle,'wb') as f:
        pickle.dump(background_dict,f,pickle.HIGHEST_PROTOCOL)
    
        
    
if __name__=="__main__":
    main()
    
