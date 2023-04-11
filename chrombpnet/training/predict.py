import pandas as pd
import os
import scipy.stats
import numpy as np
import json
import h5py
import tensorflow as tf
import chrombpnet.training.utils.argmanager as argmanager
import chrombpnet.training.utils.losses as losses
import chrombpnet.training.metrics as metrics
import chrombpnet.training.data_generators.initializers as initializers
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
from scipy import nanmean, nanstd

def write_predictions_h5py(output_prefix, profile, logcts, coords):
    # open h5 file for writing predictions
    output_h5_fname = "{}_predictions.h5".format(output_prefix)
    h5_file = h5py.File(output_h5_fname, "w")
    # create groups
    coord_group = h5_file.create_group("coords")
    pred_group = h5_file.create_group("predictions")

    num_examples=len(coords)

    coords_chrom_dset =  [str(coords[i][0]) for i in range(num_examples)]
    coords_center_dset =  [int(coords[i][1]) for i in range(num_examples)]
    coords_peak_dset =  [int(coords[i][3]) for i in range(num_examples)]

    dt = h5py.special_dtype(vlen=str)

    # create the "coords" group datasets
    coords_chrom_dset = coord_group.create_dataset(
        "coords_chrom", data=np.array(coords_chrom_dset, dtype=dt),
        dtype=dt, compression="gzip")
    coords_start_dset = coord_group.create_dataset(
        "coords_center", data=coords_center_dset, dtype=int, compression="gzip")
    coords_end_dset = coord_group.create_dataset(
        "coords_peak", data=coords_peak_dset, dtype=int, compression="gzip")

    # create the "predictions" group datasets
    profs_dset = pred_group.create_dataset(
        "profs",
        data=profile,
        dtype=float, compression="gzip")
    logcounts_dset = pred_group.create_dataset(
        "logcounts", data=logcts,
        dtype=float, compression="gzip")

    # close hdf5 file
    h5_file.close()


def load_model_wrapper(args):
    # read .h5 model
    custom_objects={"tf": tf, "multinomial_nll":losses.multinomial_nll}    
    get_custom_objects().update(custom_objects)    
    model=load_model(args.model_h5, compile=False)
    print("got the model")
    #model.summary()
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
        true_counts.extend(y[0])
        profile_probs_predictions.extend(softmax(preds[0]))

        # get profile predictions
        true_counts_sum.extend(y[1][:,0])
        counts_sum_predictions.extend(preds[1][:,0])
        coordinates.extend(coords)

    return np.array(true_counts), np.array(profile_probs_predictions), np.array(true_counts_sum), np.array(counts_sum_predictions), np.array(coordinates)


def main(args):


    metrics_dictionary = {"counts_metrics":{}, "profile_metrics":{}}

    # get model architecture to load - can load .hdf5 and .weights/.arch
    model=load_model_wrapper(args)


    test_generator = initializers.initialize_generators(args, mode="test", parameters=None, return_coords=True)
    true_counts, profile_probs_predictions, true_counts_sum, counts_sum_predictions, coordinates = predict_on_batch_wrapper(model, test_generator)


    # generate prediction on test set and store metrics
    write_predictions_h5py(args.output_prefix, profile_probs_predictions, counts_sum_predictions, coordinates)

    # store regions, their predictions and corresponding pointwise metrics
    mnll_pw, mnll_norm, jsd_pw, jsd_norm, jsd_rnd, jsd_rnd_norm, mnll_rnd, mnll_rnd_norm =  metrics.profile_metrics(true_counts,profile_probs_predictions)

    # including both metrics    
    if args.peaks != "None" and args.nonpeaks != "None":
        spearman_cor, pearson_cor, mse = metrics.counts_metrics(true_counts_sum, counts_sum_predictions,args.output_prefix+"_peaks_and_nonpeaks", "Both peaks and non peaks")
        metrics_dictionary["counts_metrics"]["peaks_and_nonpeaks"] = {}
        metrics_dictionary["counts_metrics"]["peaks_and_nonpeaks"]["spearmanr"] = spearman_cor
        metrics_dictionary["counts_metrics"]["peaks_and_nonpeaks"]["pearsonr"] = pearson_cor
        metrics_dictionary["counts_metrics"]["peaks_and_nonpeaks"]["mse"] = mse

        metrics_dictionary["profile_metrics"]["peaks_and_nonpeaks"] = {}
        metrics_dictionary["profile_metrics"]["peaks_and_nonpeaks"]["median_jsd"] = np.nanmedian(jsd_pw)        
        metrics_dictionary["profile_metrics"]["peaks_and_nonpeaks"]["median_norm_jsd"] = np.nanmedian(jsd_norm)

        metrics.plot_histogram(jsd_pw, jsd_rnd, args.output_prefix+"_peaks_and_nonpeaks", "Both peaks and non peaks")


    # including only nonpeak metrics
    if args.nonpeaks != "None":
        non_peaks_idx = coordinates[:,3] == '0'
        spearman_cor, pearson_cor, mse = metrics.counts_metrics(true_counts_sum[non_peaks_idx], counts_sum_predictions[non_peaks_idx],args.output_prefix+"_only_nonpeaks", "Only non peaks")
        metrics_dictionary["counts_metrics"]["nonpeaks"] = {}
        metrics_dictionary["counts_metrics"]["nonpeaks"]["spearmanr"] = spearman_cor
        metrics_dictionary["counts_metrics"]["nonpeaks"]["pearsonr"] = pearson_cor
        metrics_dictionary["counts_metrics"]["nonpeaks"]["mse"] = mse

        metrics_dictionary["profile_metrics"]["nonpeaks"] = {}
        metrics_dictionary["profile_metrics"]["nonpeaks"]["median_jsd"] = np.nanmedian(jsd_pw[non_peaks_idx])        
        metrics_dictionary["profile_metrics"]["nonpeaks"]["median_norm_jsd"] = np.nanmedian(jsd_norm[non_peaks_idx])

        metrics.plot_histogram(jsd_pw[non_peaks_idx], jsd_rnd[non_peaks_idx], args.output_prefix+"_only_nonpeaks", "Only non peaks")

    # including only peak metrics
    if args.peaks != "None":
        peaks_idx = coordinates[:,3] == '1'
        spearman_cor, pearson_cor, mse = metrics.counts_metrics(true_counts_sum[peaks_idx], counts_sum_predictions[peaks_idx],args.output_prefix+"_only_peaks", "Only peaks")
        metrics_dictionary["counts_metrics"]["peaks"] = {}
        metrics_dictionary["counts_metrics"]["peaks"]["spearmanr"] = spearman_cor
        metrics_dictionary["counts_metrics"]["peaks"]["pearsonr"] = pearson_cor
        metrics_dictionary["counts_metrics"]["peaks"]["mse"] = mse

        metrics_dictionary["profile_metrics"]["peaks"] = {}
        metrics_dictionary["profile_metrics"]["peaks"]["median_jsd"] = np.nanmedian(jsd_pw[peaks_idx])        
        metrics_dictionary["profile_metrics"]["peaks"]["median_norm_jsd"] = np.nanmedian(jsd_norm[peaks_idx])
        metrics.plot_histogram(jsd_pw[peaks_idx], jsd_rnd[peaks_idx], args.output_prefix+"_only_peaks", "Only peaks")

        #ofile = open(args.output_prefix+"_pearson_cor.txt","w")
        #ofile.write(str(round(pearson_cor,2)))
        #ofile.close()

        #ofile = open(args.output_prefix+"_norm_jsd.txt","w")
        #ofile.write(str(round(metrics_dictionary["profile_metrics"]["peaks"]["median_norm_jsd"],2)))
        #ofile.close()
    # store dictionary
    with open(args.output_prefix+'_metrics.json', 'w') as fp:
            json.dump(metrics_dictionary, fp,  indent=4)

if __name__=="__main__":
    # read arguments
    args=argmanager.fetch_predict_args()
    main(args)

