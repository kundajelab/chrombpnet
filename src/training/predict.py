import utils.argmanager as argmanager
import utils.losses as losses
from keras.utils.generic_utils import get_custom_objects
from tensorflow.keras.models import load_model
from scipy import nanmean, nanstd
import pandas as pd
import splits
import os
import scipy.stats
import metrics
import data_generators.initializers as initializers
import numpy as np
import json
import h5py
import tensorflow as tf

def load_model_wrapper(args):
    # read .h5 model
    custom_objects={"MultichannelMultinomialNLL": losses.MultichannelMultinomialNLL, "tf": tf}    
    get_custom_objects().update(custom_objects)    
    model=load_model(args.model_h5)
    print("got the model")
    model.summary()
    return model

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def predict_on_batch_wrapper(model,test_generator):
    num_batches=len(test_generator)
    profile_probs_predictions = []
    true_counts = []
    counts_sum_predictions = []
    true_counts_sum = []
    coordinates = []

    for idx in range(num_batches):
        if idx%100==0:
            print(str(idx)+'/'+str(num_batches))
        
        X,y,coords=test_generator[idx]

        #get the model predictions            
        preds=model.predict_on_batch(X)

        # get counts predictions
        true_counts.extend(np.squeeze(y[0]))
        profile_probs_predictions.extend(softmax(np.squeeze(preds[0])))

        # get profile predictions
        true_counts_sum.extend(np.squeeze(y[1]))
        counts_sum_predictions.extend(np.squeeze(preds[1]))
        coordinates.extend(coords)

    return np.array(true_counts), np.array(profile_probs_predictions), np.array(true_counts_sum), np.array(counts_sum_predictions), np.array(coordinates)

def get_model_param_dict(param_file):
    '''
    param_file has 2 columns -- param name in column 1, and param value in column 2
    You can pass model specfic parameters to design your own model with this
    '''
    params={}
    if param_file is None:
        return  params
    for line in open(param_file,'r').read().strip().split('\n'):
        tokens=line.split('\t')
        params[tokens[0]]=tokens[1]
    return params 

def main():

    metrics_dictionary = {"counts_metrics":{}, "profile_metrics":{}}
    # read arguments
    args=argmanager.fetch_predict_args()

    # get model architecture to load - can load .hdf5 and .weights/.arch
    model=load_model_wrapper(args)

    #parameters = get_model_param_dict(args.params)
    parameters=None

    test_generator = initializers.initialize_generators(args, mode="test", parameters=parameters, return_coords=True)

    # generate prediction on test set and store metrics
    true_counts, profile_probs_predictions, true_counts_sum, counts_sum_predictions, coordinates = predict_on_batch_wrapper(model, test_generator)

    #true_counts = scipy.ndimage.gaussian_filter1d(true_counts, 7,axis=1, truncate=(80 / 14))
 
    # store regions, their predictions and corresponding pointwise metrics
    mnll_pw, mnll_norm, jsd_pw, jsd_norm, jsd_rnd, jsd_rnd_norm, mnll_rnd, mnll_rnd_norm =  metrics.profile_metrics(true_counts,profile_probs_predictions)

    # including both metrics    
    if args.peaks != "None" and args.nonpeaks != "None":
        spearman_cor, pearson_cor, mse = metrics.counts_metrics(true_counts_sum, counts_sum_predictions,args.output_prefix+"/peaks_and_nonpeaks", "Both peaks and non peaks")
        metrics_dictionary["counts_metrics"]["spearmanr_peaks_and_nonpeaks"] = spearman_cor
        metrics_dictionary["counts_metrics"]["pearsonr_peaks_and_nonpeaks"] = pearson_cor
        metrics_dictionary["counts_metrics"]["mse_peaks_and_nonpeaks"] = mse

        metrics_dictionary["profile_metrics"]["median_jsd_peaks_and_nonpeaks"] = np.nanmedian(jsd_pw)        
        metrics_dictionary["profile_metrics"]["median_normjsd_peaks_and_nonpeaks"] = np.nanmedian(jsd_norm)
        #metrics_dictionary["profile_metrics"]["median_mnll_peaks_and_nonpeaks"] = np.median(mnll_pw)
        #metrics_dictionary["profile_metrics"]["median_normmnll_peaks_and_nonpeaks"] = np.median(mnll_norm)
   
        #metrics.plot_histogram(mnll_pw, mnll_rnd, jsd_pw, jsd_rnd, args.output_prefix+"/peaks_and_nonpeaks", "Both peaks and non peaks")
        metrics.plot_histogram(jsd_pw, jsd_rnd, args.output_prefix+"/peaks_and_nonpeaks", "Both peaks and non peaks")


    # including only nonpeak metrics
    if args.nonpeaks != "None":
        non_peaks_idx = coordinates[:,3] == '0'
        spearman_cor, pearson_cor, mse = metrics.counts_metrics(true_counts_sum[non_peaks_idx], counts_sum_predictions[non_peaks_idx],args.output_prefix+"/only_nonpeaks", "Only non peaks")
        metrics_dictionary["counts_metrics"]["spearmanr_nonpeaks"] = spearman_cor
        metrics_dictionary["counts_metrics"]["pearsonr_nonpeaks"] = pearson_cor
        metrics_dictionary["counts_metrics"]["mse_nonpeaks"] = mse

        metrics_dictionary["profile_metrics"]["median_jsd_nonpeaks"] = np.nanmedian(jsd_pw[non_peaks_idx])        
        metrics_dictionary["profile_metrics"]["median_normjsd_nonpeaks"] = np.nanmedian(jsd_norm[non_peaks_idx])

        #metrics_dictionary["profile_metrics"]["median_random_jsd_nonpeaks"] = np.nanmedian(jsd_rnd[non_peaks_idx])        
        #metrics_dictionary["profile_metrics"]["median_random_normjsd_nonpeaks"] = np.nanmedian(jsd_rnd_norm[non_peaks_idx])

        #metrics_dictionary["profile_metrics"]["median_mnll_nonpeaks"] = np.median(mnll_pw[non_peaks_idx])
        #metrics_dictionary["profile_metrics"]["median_normmnll_nonpeaks"] = np.median(mnll_pw[non_peaks_idx])

        #metrics.plot_histogram(mnll_pw[non_peaks_idx], mnll_rnd[non_peaks_idx], jsd_pw[non_peaks_idx], jsd_rnd[non_peaks_idx], args.output_prefix+"/only_nonpeaks", "Only non peaks")
        metrics.plot_histogram(jsd_pw[non_peaks_idx], jsd_rnd[non_peaks_idx], args.output_prefix+"/only_nonpeaks", "Only non peaks")

    # including only peak metrics
    if args.peaks != "None":
        peaks_idx = coordinates[:,3] == '1'
        spearman_cor, pearson_cor, mse = metrics.counts_metrics(true_counts_sum[peaks_idx], counts_sum_predictions[peaks_idx],args.output_prefix+"/only_peaks", "Only peaks")
        metrics_dictionary["counts_metrics"]["spearmanr_peaks"] = spearman_cor
        metrics_dictionary["counts_metrics"]["pearsonr_peaks"] = pearson_cor
        metrics_dictionary["counts_metrics"]["mse_peaks"] = mse

        metrics_dictionary["profile_metrics"]["median_jsd_peaks"] = np.nanmedian(jsd_pw[peaks_idx])        
        metrics_dictionary["profile_metrics"]["median_normjsd_peaks"] = np.nanmedian(jsd_norm[peaks_idx])

        #metrics_dictionary["profile_metrics"]["median_mnll_peaks"] = np.median(mnll_pw[peaks_idx])
        #metrics_dictionary["profile_metrics"]["median_normmnll_peaks"] = np.median(mnll_pw[peaks_idx])

        #metrics.plot_histogram(mnll_pw[peaks_idx], mnll_rnd[peaks_idx], jsd_pw[peaks_idx], jsd_rnd[peaks_idx], args.output_prefix+"/only_peaks", "Only peaks")
        metrics.plot_histogram(jsd_pw[peaks_idx], jsd_rnd[peaks_idx], args.output_prefix+"/only_peaks", "Only peaks")

    # store dictionary
    with open(args.output_prefix+'/metrics.json', 'w') as fp:
            json.dump(metrics_dictionary, fp,  indent=4)

    # store region-wise values
    outf=h5py.File(args.output_prefix+"region_wise_predictions_and_metrics.h5",'w')
    outf.create_dataset("true_counts",data=true_counts)
    outf.create_dataset("profile_probs_predictions",data=profile_probs_predictions)
    outf.create_dataset("true_counts_sum",data=true_counts_sum)
    outf.create_dataset("counts_sum_predictions",data=counts_sum_predictions)
    #outf.create_dataset("coordinates",data=coordinates)
    outf.close()

if __name__=="__main__":
    main()
