import argparse
import tensorflow
from tensorflow.compat.v1.keras.backend import get_session
tensorflow.compat.v1.disable_v2_behavior()
import math
import pysam
import random 
import pandas as pd
from scipy.special import softmax,expit
import pickle

#import the kerasAC dependencies
#load the model!
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import *
from kerasAC.custom_losses import *
from kerasAC.interpret.ism.ism_profile import * 
from kerasAC.helpers.transform_bpnet_io import * 
from kerasAC.util import *
from kerasAC.get_model import * 

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument('--model')
    parser.add_argument('--peak_start',type=int,default=None)
    parser.add_argument('--peak_end',type=int,default=None)
    parser.add_argument('--narrowPeak')
    parser.add_argument('--ref_fasta')
    parser.add_argument('--out_prefix')
    parser.add_argument("--flank",type=int,default=673)
    parser.add_argument("--n_sample_for_background",type=int,default=30000)
    parser.add_argument("--dump_chunk_interval",type=int,default=250) 
    return parser.parse_args() 

def get_ism_single_bp(model,seq,prof_pred,count_pred):
    #expand default preds to match dimensions for ISM 
    default_prof_expanded=np.zeros((len(seq),1000,4))
    default_count_expanded=np.zeros((len(seq),4))
    for j in range(4): 
        for i in range(len(seq)): 
            default_prof_expanded[i,:,j]=prof_pred 
            default_count_expanded[i,j]=count_pred            
    #create placeholders for ISM predictions        
    ind_to_base={0:'A',1:'C',2:'G',3:'T'}
    placeholder_prof=np.zeros((len(seq),1000,4))
    placeholder_count=np.zeros((len(seq),4))
    for j in range(4):
        cur_allele_seqs=[]
        for i in range(len(seq)): 
            cur_allele_seqs.append(seq[0:i]+ind_to_base[j]+seq[i+1::])
        #get predictions for this allele 
        cur_allele_preds=model.predict([one_hot_encode(cur_allele_seqs)])
        cur_allele_prof=np.squeeze(cur_allele_preds[0])
        cur_allele_count=np.squeeze(cur_allele_preds[1])
        placeholder_prof[:,:,j]=cur_allele_prof
        placeholder_count[:,j]=cur_allele_count
    
    #subtract the WT, average across base axis 
    placeholder_prof_normed=placeholder_prof-default_prof_expanded
    placeholder_count_normed=placeholder_count-default_count_expanded

    placeholder_prof_normed=placeholder_prof_normed-np.expand_dims(np.mean(placeholder_prof_normed,axis=2),axis=2)
    placeholder_count_normed=placeholder_count_normed-np.expand_dims(np.mean(placeholder_count_normed,axis=1),axis=1)
    
    seq_onehot=one_hot_encode([seq])
    
    ism_profile_track=np.sum(np.abs(placeholder_prof_normed),axis=1)
    ism_count_track=placeholder_count_normed
    
    #ism heatmap: 
    ism_mat_observed=placeholder_prof_normed
    return ism_profile_track, ism_count_track, ism_mat_observed

def get_model(hdf5_file):
    #load the model!
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
    model=load_model(hdf5_file)
    return model


def main():
    args=parse_args()
    #load the model
    model=get_model(args.model) 
    print("loaded the model")
    peaks=pd.read_csv(args.narrowPeak,header=None,sep="\t")
    print("loaded narrowPeak")
    
    #get the reference fasta file
    ref=pysam.FastaFile(args.ref_fasta)

    #create placeholders for background tracks
    background_profile_ism=np.empty([args.flank*2*args.n_sample_for_background])
    background_counts_ism=np.empty([args.flank*2*args.n_sample_for_background])
    background_index_interval=peaks.shape[0]/args.n_sample_for_background

    ism_foreground_score_dict={}
    if args.peak_start is None:
        peak_start=0
    else:
        peak_start=args.peak_start
    if args.peak_end is None:
        peak_end=peaks.shape[0]
    else:
        peak_end=args.peak_end
        
    for index,row in peaks.iterrows():
        if index < peak_start:
            continue
        if index > peak_end:
            continue
        print(str(index))
        chrom=row[0]
        summit=row[1]+row[9]
        seq=ref.fetch(chrom,summit-args.flank,summit+args.flank)
        preds_seq=model.predict([one_hot_encode([seq])])
        preds_seq_prof=np.squeeze(preds_seq[0])
        preds_seq_count=np.squeeze(preds_seq[1])
        #get ISM scores 
        seq_ism_profile_track, seq_ism_count_track, seq_ism_mat_observed = get_ism_single_bp(model,seq,preds_seq_prof,preds_seq_count)
        ism_foreground_score_dict[(chrom,summit,row[1],row[2])]=[seq_ism_profile_track,seq_ism_count_track,seq_ism_mat_observed]
        
        if index % background_index_interval==0: 
            #get predictions for the sequence and the scrambled sequence
            scrambled_seq=''.join(random.sample(seq,len(seq))) #generate scrambled sequence to calculate background ISM score distribution
            scrambled_one_hot=one_hot_encode([scrambled_seq])
            preds_scrambled_seq=model.predict([scrambled_one_hot])
            preds_scrambled_seq_prof=np.squeeze(preds_scrambled_seq[0])
            preds_scrambled_seq_count=np.squeeze(preds_scrambled_seq[1])
            #get ISM scores for scrambled track
            scrambled_ism_profile_track, scrambled_ism_count_track,scrambled_ism_mat_observed = get_ism_single_bp(model,scrambled_seq,preds_scrambled_seq_prof,preds_scrambled_seq_count)
            #flatten the scores for background
            background_profile_ism[index*(2*args.flank):(index+1)*(2*args.flank)]=np.sum(np.squeeze(scrambled_one_hot*scrambled_ism_profile_track),axis=1)
            background_counts_ism[index*(2*args.flank):(index+1)*(2*args.flank)]=np.sum(np.squeeze(scrambled_one_hot*scrambled_ism_count_track),axis=1)
            
    #save to output files
    pickle.dump( background_profile_ism, open( args.out_prefix+"background_profile_ism.p", "wb" ) )
    print("dumped background profile ISM to pickle")
    pickle.dump(background_counts_ism,open(args.out_prefix+"background_counts_ism.p","wb"))
    print("dumpted background counts ISM to pickle")
    pickle.dump(ism_foreground_score_dict,open(args.out_prefix+"ism_foreground_score_dict.p","wb"))
        
        
if __name__=="__main__":
    main()
    
