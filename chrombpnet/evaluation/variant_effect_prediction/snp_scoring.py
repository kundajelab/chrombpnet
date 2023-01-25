from chrombpnet.evaluation.variant_effect_prediction.snp_generator import SNPGenerator
from scipy.spatial.distance import jensenshannon
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import chrombpnet.training.utils.losses as losses
import pandas as pd
import os
import argparse
import numpy as np
import pickle as pkl
import tensorflow as tf


SNP_SCHEMA = ["CHR", "POS0", "REF", "ALT", "META_DATA"]

def fetch_variant_args():
    parser=argparse.ArgumentParser(description="variant effect scoring scripts on SNPS")
    parser.add_argument("-snps", "--snp_data", type=str, required=True, help="Path to a tsv output with the following information in columns - chr, position to insert allele (0-based), ref allele, alt allele")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-m","--model_h5", type=str, required=True, help="Path to model hdf5")
    parser.add_argument("-op","--output-prefix", type=str, required=True, help="Path to storing snp effect score predictions from the script, directory should already exist")
    parser.add_argument("-bs","--batch_size", type=int, default=64, help="Batch size to use for model")
    parser.add_argument("-dm","--debug_mode_on", type=int, default=0, help="Use this mode to print the flanks of first five SNP insert locations")
    args = parser.parse_args()
    return args

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def load_model_wrapper(args):
    # read .h5 model
    custom_objects={"tf": tf, "multinomial_nll":losses.multinomial_nll}    
    get_custom_objects().update(custom_objects)    
    model=load_model(args.model_h5)
    print("model loaded succesfully")
    return model

def fetch_snp_predictions(model, snp_regions, inputlen, genome_fasta, batch_size, debug_mode_on=False):
    '''
    Returns model predictions (counts and profile probability predictions) at the given reference and alternate snp alleles.
    Please note that if the SNP location is at the edge - i.e we are unable to form a given inputlen of sequence - we skip predictions at this SNP

    Arguments::
        model: chrombpnet model .h5 file to use for snp scoring. 
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
    rsids = []
    ref_logcount_preds=[]
    alt_logcount_preds=[]
    ref_prob_preds=[]
    alt_prob_preds=[]

    # snp sequence generator 
    snp_gen=SNPGenerator(snp_regions=snp_regions,
                        inputlen=inputlen,
                        genome_fasta=genome_fasta,
                        batch_size=batch_size,
                        debug_mode_on=debug_mode_on)

    for i in range(len(snp_gen)):

        batch_rsids, ref_seqs, alt_seqs = snp_gen[i]

        ref_batch_preds=model.predict(ref_seqs)
        alt_batch_preds=model.predict(alt_seqs)

        ref_logcount_preds.extend(np.squeeze(ref_batch_preds[1]))
        alt_logcount_preds.extend(np.squeeze(alt_batch_preds[1]))

        ref_prob_preds.extend(np.squeeze(softmax(ref_batch_preds[0])))
        alt_prob_preds.extend(np.squeeze(softmax(alt_batch_preds[0])))

        rsids.extend(batch_rsids)

    return np.array(rsids), np.array(ref_logcount_preds), np.array(alt_logcount_preds), np.array(ref_prob_preds), np.array(alt_prob_preds)

def predict_snp_effect_scores(rsids, ref_count_preds, alt_count_preds, ref_prob_preds, alt_prob_preds):
    '''
    Predicts variant effect scores based on model predictions.

    Arguments::
       ref_logcount_preds: log count predictions at the reference allele with size (N,)
       alt_logcount_preds: log count predictions at the alternate alele with size (N,)
       ref_prob_preds: profile probability predictions at the reference allele with size (N,outputlen). outputlen depends on the model.
       alt_prob_preds:  profile probability predictions at the alternate allele with size (N,outputlen). outputlen depends on the model.
    
    Returns:
        log_counts_diff: difference in log count predictions of alternate and reference allele (N,)
        log_probs_diff_abs_sum: Sum of absolute difference in log probability prediction of alternate and reference allele per base. (N,)
        probs_jsd_diff: Jensenshannon distance between probability predictions of alternate and reference allele (N,)
    '''
    log_counts_diff = alt_count_preds - ref_count_preds
    log_probs_diff_abs_sum =  np.sum(np.abs(np.log(alt_prob_preds) -  np.log(ref_prob_preds)),axis=1)*np.sign(log_counts_diff)
    probs_jsd_diff = np.array([jensenshannon(x,y) for x,y in zip(alt_prob_preds, ref_prob_preds)])*np.sign(log_counts_diff)

    return log_counts_diff, log_probs_diff_abs_sum, probs_jsd_diff


def main(args):
    debug_mode_on = args.debug_mode_on

    # load the model
    model = load_model_wrapper(args)

    # load the snp data
    snp_regions=pd.read_csv(args.snp_data,header=None,sep='\t', names=SNP_SCHEMA)
    snp_regions["META_DATA"].fillna('', inplace=True)
    snp_regions['RSID']=snp_regions['CHR'].astype(str)+'_'+snp_regions['POS0'].astype(str)+'_'+snp_regions['REF'].astype(str)+'_'+snp_regions['ALT'].astype('str')+"_"+snp_regions['META_DATA'].astype('str')
    print("printing first 5 rows of the input SNP data provided..")
    print(snp_regions.head(5))

    if debug_mode_on:
        snp_regions = snp_regions.head(5)

    # infer input length
    inputlen=model.input_shape[1]
    print("input length inferred from the model: ", inputlen)

    # fetch model prediction on snps
    rsids, ref_logcount_preds, alt_logcount_preds, ref_prob_preds, alt_prob_preds = fetch_snp_predictions(model, snp_regions, inputlen, args.genome, args.batch_size, debug_mode_on)

    # find varaint effect scores at snps
    log_counts_diff, log_probs_diff_abs_sum, probs_jsd_diff = predict_snp_effect_scores(rsids, ref_logcount_preds, alt_logcount_preds, ref_prob_preds, alt_prob_preds)

    # unpack rsids to write outputs and write score to output
    snp_effect_scores_pd=pd.DataFrame()
    snp_effect_scores_pd[["CHR", "POS0", "REF", "ALT", "META_DATA"]] = pd.Series(rsids).str.split('_', expand=True)
    snp_effect_scores_pd["log_counts_diff"] = log_counts_diff
    snp_effect_scores_pd["log_probs_diff_abs_sum"] = log_probs_diff_abs_sum
    snp_effect_scores_pd["probs_jsd_diff"] = probs_jsd_diff

    snp_effect_scores_pd.to_csv(args.output_prefix+"_snp_scores.tsv", sep="\t", index=False)

    # store predictions at snps too - can compute variant effect metrics of your interest - let me know if you find something interesting :)
    data={}
    data["rsids"] = rsids
    data["ref_logcount_preds"] = ref_logcount_preds
    data["alt_logcount_preds"] = alt_logcount_preds
    data["ref_prob_preds"] = ref_prob_preds
    data["alt_prob_preds"] = alt_prob_preds

    
    pkl.dump(data, open(args.output_prefix+"_predictions_at_snp.pkl",'wb'))



if __name__=="__main__":
    args = fetch_variant_args()
    main(args)

