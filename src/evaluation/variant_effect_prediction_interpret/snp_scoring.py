from snp_generator import SNPGenerator
from scipy.spatial.distance import jensenshannon
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import pandas as pd
import os
import argparse
import losses
import numpy as np
import pickle as pkl
import deepdish as dd
import json
import numpy as np
import tensorflow as tf
import pandas as pd
import shap
from tensorflow import keras
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import pyfaidx
import shutil
import errno
import os
import argparse
import shap_utils 
import context
import pickle

tf.compat.v1.disable_eager_execution()


SNP_SCHEMA = ["CHR", "POS0", "REF", "ALT", "META_DATA"]

def fetch_variant_args():
    parser=argparse.ArgumentParser(description="variant effect scoring scripts on SNPS")
    parser.add_argument("-i", "--snp_data", type=str, required=True, help="Path to a tsv output with the following information in columns - chr, position to insert allele (0-based), ref allele, alt allele")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-m","--model_h5", type=str, required=True, help="Path to model hdf5")
    parser.add_argument("-o","--output_dir", type=str, required=True, help="Path to storing snp effect score predictions from the script, directory should already exist")
    parser.add_argument("-bs","--batch_size", type=int, default=64, help="Batch size to use for model")
    parser.add_argument("-dm","--debug_mode_on", type=int, default=0, help="Use this mode to print the flanks of first five SNP insert locations")
    args = parser.parse_args()
    return args

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def load_model_wrapper(args):
    # read .h5 model
    custom_objects={"MultichannelMultinomialNLL": losses.MultichannelMultinomialNLL}    
    get_custom_objects().update(custom_objects)    
    model=load_model(args.model_h5)
    print("model loaded succesfully")
    return model


def interpret(model, seqs, profile_or_counts):
    print("Seqs dimension : {}".format(seqs.shape))

    outlen = model.output_shape[0][1]

    profile_model_input = model.input
    profile_input = seqs
    counts_model_input = model.input
    counts_input = seqs

    if "counts" in profile_or_counts:
        profile_model_counts_explainer = shap.explainers.deep.TFDeepExplainer(
            (counts_model_input, tf.reduce_sum(model.outputs[1], axis=-1)),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref)

        print("Generating 'counts' shap scores")
        counts_shap_scores = profile_model_counts_explainer.shap_values(
            counts_input, progress_message=100)
    
        return counts_shap_scores

    if "profile" in profile_or_counts:
        weightedsum_meannormed_logits = shap_utils.get_weightedsum_meannormed_logits(model)
        profile_model_profile_explainer = shap.explainers.deep.TFDeepExplainer(
            (profile_model_input, weightedsum_meannormed_logits),
            shap_utils.shuffle_several_times,
            combine_mult_and_diffref=shap_utils.combine_mult_and_diffref)

        print("Generating 'profile' shap scores")
        profile_shap_scores = profile_model_profile_explainer.shap_values(
            profile_input, progress_message=100)

        return profile_shap_scores


def fetch_snp_predictions(snp_regions, inputlen, genome_fasta, batch_size, output_path, debug_mode_on=False):
    '''
    Returns model predictions (counts and profile probability predictions) at the given reference and alternate snp alleles.
    Please note that if the SNP location is at the edge - i.e we are unable to form a given inputlen of sequence - we skip predictions at this SNP

    Arguments::
        snp_regions: pandas dataframe with the following columns "CHR", "POS0", "REF", "ALT"
        inputlen: integer representing the input length to use, snp is inserted in the middle
        genome_fasta: path to reference genome
        batch_size: integer value with batch size to use for the model
        debug_mode_on: Takes 0/1 value. Set this  to 1 to print the flanks of first five SNP insert locations. Predictions will be provided only on the these 5 locations.
    
    Returns:
       rsids: Numpy array with (N,) SNP ids. SNP id is a string with the following values "CHR", "POS0", "REF", "ALT" concatenated with delimiter "_". 
            For each of these ids we return the predictions in the lists below. 
       ref_logcount_preds: log count predictions at the reference allele with size (N,)
       alt_logcount_preds: log count predictions at the alternate alele with size (N,)
       ref_prob_preds: profile probability predictions at the reference allele with size (N,outputlen). outputlen depends on the model.
       alt_prob_preds:  profile probability predictions at the alternate allele with size (N,outputlen). outputlen depends on the model.
    '''

    # snp sequence generator 
    snp_gen=SNPGenerator(snp_regions=snp_regions,
                        inputlen=inputlen,
                        genome_fasta=genome_fasta,
                        batch_size=batch_size,
                        debug_mode_on=debug_mode_on)

    count_preds={} 
    profile_preds={} 
    profile_shap={}
    count_shap={}
    snp_to_seq={}

    for i in range(len(snp_gen)):

        batch_rsids, ref_seqs, alt_seqs = snp_gen[i]

        ref_batch_preds=model.predict(ref_seqs)
        alt_batch_preds=model.predict(alt_seqs)

        ref_profile_interpret = interpret(model, ref_seqs,  ["profile"])
        ref_counts_interpret = interpret(model, ref_seqs,  ["counts"])

        alt_profile_interpret = interpret(model, alt_seqs,  ["profile"])
        alt_counts_interpret = interpret(model, alt_seqs,  ["counts"])


        batch_preds_profile=ref_batch_preds[0]
        batch_preds_count=ref_batch_preds[1] 
        seqs=ref_seqs
        profile_explanations=ref_profile_interpret
        count_explanations=ref_counts_interpret

        for batch_index in range(len(batch_rsids)): 
            cur_rsid=batch_rsids[batch_index]
    
            cur_pred_profile=batch_preds_profile[batch_index,:]
            cur_pred_count=batch_preds_count[batch_index,:]
            count_preds[cur_rsid]={}
            count_preds[cur_rsid]['ref']=cur_pred_count[0]
            profile_preds[cur_rsid]={}
            profile_preds[cur_rsid]['ref']=cur_pred_profile 
            profile_shap[cur_rsid]={}
            profile_shap[cur_rsid]['ref']=np.sum(profile_explanations[batch_index,:,:]*seqs[batch_index,:,:], axis=1)
            count_shap[cur_rsid]={}
            count_shap[cur_rsid]['ref']=np.sum(count_explanations[batch_index,:,:]*seqs[batch_index,:,:], axis=1)
            print("ref profile")
            print(np.max(profile_shap[cur_rsid]['ref']))
            print("ref counts")
            print(np.max(count_shap[cur_rsid]['ref']))


        batch_preds_profile=alt_batch_preds[0]
        batch_preds_count=alt_batch_preds[1] 
        seqs=alt_seqs
        profile_explanations=alt_profile_interpret
        count_explanations=alt_counts_interpret

        for batch_index in range(len(batch_rsids)): 
            cur_rsid=batch_rsids[batch_index]
    
            cur_pred_profile=batch_preds_profile[batch_index,:]
            cur_pred_count=batch_preds_count[batch_index,:]
            count_preds[cur_rsid]['alt']=cur_pred_count[0]
            profile_preds[cur_rsid]['alt']=cur_pred_profile 
            profile_shap[cur_rsid]['alt']=np.sum(profile_explanations[batch_index,:,:]*seqs[batch_index,:,:], axis=1)
            count_shap[cur_rsid]['alt']=np.sum(count_explanations[batch_index,:,:]*seqs[batch_index,:,:], axis=1)
            print("alt profile")
            print(np.max(profile_shap[cur_rsid]['alt']))
            print("alt count")
            print(np.max(count_shap[cur_rsid]['alt']))



    #Store counts
    count_preds_df=pd.DataFrame.from_dict(count_preds,orient='index')
    count_preds_df["rsid"] = count_preds_df.index.values
    count_preds_df.to_csv(os.path.join(output_path, "count_predictions_alt_and_ref.tsv"),header=True,index=True,sep='\t')

    #Store profile preds 
    pickle.dump(profile_preds, open(os.path.join(output_path,"profile_predictions.pkl"), "wb" )) 

    #Store profile shap 
    pickle.dump(profile_shap, open(os.path.join(output_path,"profile_shap.pkl"), "wb" )) 

    #Store count shap
    pickle.dump(count_shap, open(os.path.join(output_path,"count_shap.pkl"), "wb" )) 


if __name__=="__main__":

    args = fetch_variant_args()
    debug_mode_on = args.debug_mode_on

    # load the model
    model = load_model_wrapper(args)

    # load the snp data
    snp_regions=pd.read_csv(args.snp_data,header=None,sep='\t', names=SNP_SCHEMA)
    snp_regions["META_DATA"].fillna('', inplace=True)
    snp_regions['POS0'] = snp_regions['POS0'] - 1
    snp_regions['RSID']=snp_regions['CHR'].astype(str)+'_'+snp_regions['POS0'].astype(str)+'_'+snp_regions['REF'].astype(str)+'_'+snp_regions['ALT'].astype('str')+"_"+snp_regions['META_DATA'].astype('str')
    print("printing first 5 rows of the input SNP data provided..")
    print(snp_regions.head(5))

    if debug_mode_on:
        snp_regions = snp_regions.head(5)

    # infer input length
    inputlen=model.input_shape[1]
    print("input length inferred from the model: ", inputlen)

    # fetch model prediction on snps
    fetch_snp_predictions(snp_regions, inputlen, args.genome, args.batch_size, args.output_dir, debug_mode_on)


