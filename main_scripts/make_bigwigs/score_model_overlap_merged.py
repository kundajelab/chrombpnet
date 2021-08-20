import argparse
import pickle
import operator
import pysam
import numpy as np 
from scipy.special import softmax 
import tensorflow
import pybedtools
import pyBigWig
from tensorflow.compat.v1.keras.backend import get_session
tensorflow.compat.v1.disable_v2_behavior()
import math
import kerasAC 
from scipy.special import softmax,expit
from kerasAC.interpret.deepshap import * 
from kerasAC.interpret.profile_shap import * 
from kerasAC.vis import * 
from kerasAC.helpers.transform_bpnet_io import * 
from kerasAC.util import * 
import pandas as pd
from statsmodels.distributions.empirical_distribution import ECDF
import pyBigWig
import pdb 
#load the model! 
from keras.models import load_model
from keras.utils.generic_utils import get_custom_objects
from kerasAC.metrics import * 
from kerasAC.custom_losses import *
from load_model import *

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--model_path", default=None)
    parser.add_argument("--json_string", default=None)
    parser.add_argument("--weights", default=None)
    parser.add_argument("--batch_size",type=int,default=1000)
    parser.add_argument("--bigwig_labels")
    parser.add_argument("--out_prefix")
    parser.add_argument("--bed_file_to_score")
    parser.add_argument("--precentered_intervals",action='store_true',default=False,help="bed file contains intervals to be scored directly")
    parser.add_argument("--presorted_intervals",action="store_true",default=False,help="if flag provided, indicates bed file is already sorted")
    parser.add_argument("--flank_input",type=int,default=1057)
    parser.add_argument("--flank_output",type=int,default=500) 
    parser.add_argument("--ref_fasta",default="/mnt/data/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta")
    parser.add_argument("--chrom_sizes")
    parser.add_argument("--run_ism",default=False,action="store_true")
    parser.add_argument("--batch-size",type=int,default=500)
    return parser.parse_args()

def get_deepshap(prof_explainer,count_explainer,seq_onehot):
    #dim will be (1,2114,4)
    hyp_profile_explanations=prof_explainer.shap_values(seq_onehot)
    hyp_count_explanations=count_explainer.shap_values(seq_onehot)[0]
    
    profile_explanations=np.squeeze(hyp_profile_explanations*seq_onehot).sum(axis=-1)
    count_explanations=np.squeeze(hyp_count_explanations*seq_onehot).sum(axis=-1)
    
    return np.squeeze(hyp_profile_explanations), np.squeeze(hyp_count_explanations), profile_explanations, count_explanations 


def get_preds(model,seq_onehot):
    preds=model.predict(seq_onehot)
    prof=preds[0] #logits
    probs=np.squeeze(softmax(prof,axis=1)) #probabilities 
    count=np.squeeze(preds[1])   #count (1 val) 
    count_track=probs*np.expand_dims(np.exp(count),axis=1)  #count track 
    #print(np.exp(count))
    return prof,count,probs,count_track


# 2114, 1000, 4
#get observed ism, sum by pos & neg and take the one w/ larger magnitude 
def get_observed_ism_from_mat(profile,seq):
    pos_mask=np.zeros(profile.shape)
    pos_mask[profile>0]=1
    neg_mask=np.zeros(profile.shape)
    neg_mask[profile<0]=1
    
    pos_sum=np.sum(profile*pos_mask,axis=1)
    neg_sum=np.sum(profile*neg_mask,axis=1)
    
    pos_sum_mask=np.zeros(pos_sum.shape)
    pos_sum_mask[pos_sum>abs(neg_sum)]=1
    neg_sum_mask=np.zeros(neg_sum.shape)
    neg_sum_mask[neg_sum<0]=1

    profile_summed_ism=pos_sum*pos_sum_mask+neg_sum*neg_sum_mask
    
    profile_summed_scaled=np.sum(np.squeeze(profile_summed_ism*seq),axis=-1)  
    return profile_summed_scaled

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

    return placeholder_prof_normed, placeholder_count_normed


def main():
    args=parse_args()
    if args.model_path == None:
        model=load_model_wrapper(json_string=args.json_string, weights=args.weights)
    else:        
        model=load_model_wrapper(model_hdf5=args.model_path)

    print("loaded model") 
    pbw=pyBigWig.open(args.bigwig_labels,'r')

    #create the explainers 
    model_wrapper=(model.input, model.outputs[1])
    count_explainer=shap.DeepExplainer(model_wrapper,
                                   data=create_background_atac,
                                   combine_mult_and_diffref=combine_mult_and_diffref_1d)
    print("created shap explainers") 
    
    prof_output = model.output[0]                                                                                                      
    logits = prof_output - tf.reduce_mean(prof_output, axis=1, keepdims=True)
    logits_stopgrad=tf.stop_gradient(logits)
    probs=tf.nn.softmax(logits_stopgrad,axis=1)
    logits_weighted=logits*probs 
    prof_sum=tf.reduce_sum(logits_weighted,axis=(1,2))
    prof_explainer=shap.DeepExplainer(model=(model.input,prof_sum),
                                  data=create_background_atac,
                                  combine_mult_and_diffref=combine_mult_and_diffref_1d)


    ref=pysam.FastaFile(args.ref_fasta)
    regions=pd.read_csv(args.bed_file_to_score.strip(),header=None,sep='\t')
    #make sure that the regions are sorted!!!
    if args.presorted_intervals is False:
        if regions.shape[1]>9: 
            regions=regions.sort_values(by=[0,1,2,9])
        else:
            regions=regions.sort_values(by=[0,1,2])
        print("sorted input regions")

    #initialize the bigwigs
    #prepare chromsizes for bigwig header use
    chromsizes=[i.split('\t') for i in open(args.chrom_sizes.strip(),'r').read().strip().split('\n')]
    chromsizes=[tuple([i[0],int(i[1])]) for i in chromsizes]
    chromsizes=sorted(chromsizes,key=operator.itemgetter(0))
    chroms=[i[0] for i in chromsizes]
    print(str(chroms))
    print(str(chromsizes)) 
    chromsize_dict={}
    for entry in chromsizes:
        chromsize_dict[entry[0]]=int(entry[1])

    outputs={}
    
    #labels -- we already have the input in bigwig format, but helpful to restrict just to peak regions for visualization
    outputs['observed_count_bigwig']=pyBigWig.open(args.out_prefix+".observed.counts.bw",'w')
    outputs['observed_count_bigwig'].addHeader(chromsizes)

    outputs['predicted_count_bigwig']=pyBigWig.open(args.out_prefix+".predicted.counts.bw",'w')
    outputs['predicted_count_bigwig'].addHeader(chromsizes) 

    outputs['predicted_probability_bigwig']=pyBigWig.open(args.out_prefix+".predicted.probability.bw",'w')
    outputs['predicted_probability_bigwig'].addHeader(chromsizes)

    outputs['observed_probability_bigwig']=pyBigWig.open(args.out_prefix+".observed.probability.bw",'w')
    outputs['observed_probability_bigwig'].addHeader(chromsizes)
    
    outputs['deepshap_profile']=pyBigWig.open(args.out_prefix+".deepshap.profile.bw",'w')
    outputs['deepshap_profile'].addHeader(chromsizes)
    
    outputs['deepshap_counts']=pyBigWig.open(args.out_prefix+".deepshap.counts.bw",'w')
    outputs['deepshap_counts'].addHeader(chromsizes)
    
    if args.run_ism is True:
        outputs['ism_profile']=pyBigWig.open(args.out_prefix+".ism.profile.bw",'w')
        outputs['ism_profile'].addHeader(chromsizes)
        outputs['ism_count']=pyBigWig.open(args.out_prefix+".ism.count.bw",'w')
        outputs['ism_count'].addHeader(chromsizes)
        
    for chrom in chroms:
        print(str(chrom))
        chrom_dict={}
        chrom_dict['observed_count_bigwig']=np.empty((chromsize_dict[chrom],))
        chrom_dict['predicted_count_bigwig']=np.empty((chromsize_dict[chrom],))
        chrom_dict['predicted_probability_bigwig']=np.empty((chromsize_dict[chrom],))
        chrom_dict['observed_probability_bigwig']=np.empty((chromsize_dict[chrom],))
        chrom_dict['deepshap_profile']=np.empty((chromsize_dict[chrom],))
        chrom_dict['deepshap_counts']=np.empty((chromsize_dict[chrom],))
        chrom_dict['deepshap_profile_hyp']=np.empty((chromsize_dict[chrom],4))
        chrom_dict['deepshap_counts_hyp']=np.empty((chromsize_dict[chrom],4))
        if args.run_ism is True:
            chrom_dict['ism_profile']=np.empty((chromsize_dict[chrom],))
            chrom_dict['ism_count']=np.empty((chromsize_dict[chrom],))

        cur_batch_n=0
        regions_chrom=regions[regions[0]==chrom]
        num_batches=regions_chrom.shape[0]//args.batch_size +1
        
        for batch in np.array_split(regions_chrom,num_batches):
            batch=batch.reset_index(drop=True)
            cur_batch_n+=1
            print(str(cur_batch_n))
            cur_batch_size=batch.shape[0]
            if args.precentered_intervals is False:
                summit=batch[1]+batch[9]
                start_pos_input=summit-args.flank_input
                end_pos_input=summit+args.flank_input
                start_pos_output=summit-args.flank_output
                end_pos_output=summit+args.flank_output
            else:
                start_pos_input=batch[1]
                end_pos_input=batch[2]
                center=(start_pos_input+end_pos_input)//2
                start_pos_output=center-args.flank_output
                end_pos_output=center+args.flank_output


            #get the reference and alternate one-hot-encoded sequences 
            seq=[ref.fetch(chrom,start_pos_input[i],end_pos_input[i]) for i in range(cur_batch_size)]
            onehot=one_hot_encode(seq)

            #profile head out, count head out, probability predicted track, count predicted track 
            predicted_profile_logit_head,predicted_count_head,predicted_prob_track,predicted_count_track=get_preds(model,onehot)

            #get deepSHAP scores  
            hyp_profile_explanations_shap, hyp_count_explanations_shap,profile_explanations_shap, count_explanations_shap=get_deepshap(prof_explainer, count_explainer, onehot)

            #get ISM scores
            if args.run_ism ==True:
                single_bp_ism_profile_track,single_bp_ism_count_track= get_ism_single_bp(model,
                                                                                         seq,
                                                                                         np.squeeze(predicted_profile_logit_head),
                                                                                         np.squeeze(predicted_count_head))
                single_bp_ism_profile_track_adjusted=get_observed_ism_from_mat(single_bp_ism_profile_track,
                                                                               onehot)
                single_bp_ism_count_track_adjusted=np.sum(np.squeeze(single_bp_ism_count_track*onehot),axis=-1)


            #store in dictionary 
            for batch_index in range(cur_batch_size):
                cur_chrom=chrom
                cur_start_pos_output=start_pos_output[batch_index]
                cur_end_pos_output=end_pos_output[batch_index]
                cur_start_pos_input=start_pos_input[batch_index]
                cur_end_pos_input=end_pos_input[batch_index]


                chrom_dict['observed_count_bigwig'][cur_start_pos_output:cur_end_pos_output]=np.nan_to_num(pbw.values(cur_chrom,cur_start_pos_output,cur_end_pos_output))
                chrom_dict['predicted_count_bigwig'][cur_start_pos_output:cur_end_pos_output]=np.squeeze(predicted_count_track[batch_index])
                chrom_dict['predicted_probability_bigwig'][cur_start_pos_output:cur_end_pos_output]=np.squeeze(predicted_prob_track[batch_index])
                chrom_dict['observed_probability_bigwig'][cur_start_pos_output:cur_end_pos_output]= chrom_dict['observed_count_bigwig'][cur_start_pos_output:cur_end_pos_output] /  sum(chrom_dict['observed_count_bigwig'][cur_start_pos_output:cur_end_pos_output])
                chrom_dict['deepshap_profile_hyp'][cur_start_pos_input:cur_end_pos_input]=np.squeeze(hyp_profile_explanations_shap[batch_index])
                chrom_dict['deepshap_counts_hyp'][cur_start_pos_input:cur_end_pos_input]=np.squeeze(hyp_count_explanations_shap[batch_index])
                chrom_dict['deepshap_profile'][cur_start_pos_input:cur_end_pos_input]=np.squeeze(profile_explanations_shap[batch_index])
                chrom_dict['deepshap_counts'][cur_start_pos_input:cur_end_pos_input]=np.squeeze(count_explanations_shap[batch_index])
                if args.run_ism is True:
                    chrom_dict['ism_profile']=np.squeeze(single_bp_ism_profile_track_adjusted)
                    chrom_dict['ism_count']=np.squeeze(single_bp_ism_count_track_adjusted)

        #dump hypothetical scores to pickles
        #write to bigwig
        print("writing to bigwig") 
        print(str(chrom))
        #with open(args.out_prefix+'.hyp.profile.deepshap.'+chrom+'.pkl','wb') as handle:
        #    pickle.dump(chrom_dict['deepshap_profile_hyp'],handle,protocol=pickle.HIGHEST_PROTOCOL)
        #with open(args.out_prefix+'.hyp.count.deepshap.'+chrom+'.pkl','wb') as handle:
        #    pickle.dump(chrom_dict['deepshap_counts_hyp'],handle,protocol=pickle.HIGHEST_PROTOCOL) 
        print("pickled hypothetical scores") 
        outputs['observed_count_bigwig'].addEntries(chrom,0,values=chrom_dict['observed_count_bigwig'], span=1,step=1)
        outputs['predicted_count_bigwig'].addEntries(chrom,0,values=chrom_dict['predicted_count_bigwig'],span=1,step=1)
        outputs['predicted_probability_bigwig'].addEntries(chrom,0,values=chrom_dict['predicted_probability_bigwig'],span=1,step=1)
        outputs['observed_probability_bigwig'].addEntries(chrom,0,values=chrom_dict['observed_probability_bigwig'],span=1,step=1)
        outputs['deepshap_profile'].addEntries(chrom,0,values=chrom_dict['deepshap_profile'],span=1,step=1)
        outputs['deepshap_counts'].addEntries(chrom,0,values=chrom_dict['deepshap_counts'],span=1,step=1)
        if args.run_ism is True:
            outputs['ism_profile'].addEntries(chrom,0,values=chrom_dict['ism_profile'],span=1,step=1)
            outputs['ism_count'].addEntries(chrom,0,values=chrom_dict['ism_count'],span=1,step=1)
    #close the files
    for key in outputs:
        outputs[key].close() 

if __name__=="__main__":
    main()
    
