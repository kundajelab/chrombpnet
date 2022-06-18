from __future__ import division, print_function, absolute_import
import importlib.machinery
import tensorflow.keras.callbacks as tfcallbacks 
import utils.bias_argmanager as argmanager
import utils.losses as losses
import utils.callbacks as callbacks
import data_generators.tobias_initializers as initializers
import pandas as pd
import os
import json
import numpy as np
import matplotlib.pyplot as plt
from tensorflow.keras.utils import get_custom_objects
from tensorflow.keras.models import load_model
import tensorflow as tf

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
os.environ['PYTHONHASHSEED'] = '0'
import matplotlib.pyplot as plt

def load_model_wrapper(model_h5):
    # read .h5 model
    custom_objects={"tf": tf, "multinomial_nll":losses.multinomial_nll}    
    get_custom_objects().update(custom_objects)    
    model=load_model(model_h5)
    print("got the model")
    #model.summary()
    return model

def softmax(x, temp=1):
    norm_x = x - np.mean(x,axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)



def predict_on_batch_wrapper(test_generator, model):
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
    return np.array(true_counts), np.array(true_counts_sum), np.array(coordinates), np.array(profile_probs_predictions), np.array(counts_sum_predictions)



def get_model_param_dict(args):
    '''
    param_file is a TSV file with 2 columns -- param name in column 1, and param value in column 2
    You can pass model specfic parameters to design your own model with this.
    '''
    params={}
    for line in open(args.params,'r').read().strip().split('\n'):
        tokens=line.split('\t')
        params[tokens[0]]=tokens[1]

    assert("counts_loss_weight" in params.keys()) # missing counts loss weight to use
    assert("filters" in params.keys()) # filters to use for the model not provided
    assert("n_dil_layers" in params.keys()) # n_dil_layers to use for the model not provided
    assert("inputlen" in params.keys()) # inputlen to use for the model not provided
    assert("outputlen" in params.keys()) # outputlen to use for the model not provided
    assert("negative_sampling_ratio" in params.keys()) # negative_sampling_ratio to use for the model not provided
    assert("max_jitter" in params.keys()) # max_jitter to use for the model not provided
    #assert(args.chr_fold_path==params["chr_fold_path"]) # the parameters were generated on a different folds compared to the given fold

    assert(int(params["inputlen"])%2==0)
    assert(int(params["outputlen"])%2==0)

    return params 

def main():

    # read arguments
    args=argmanager.fetch_train_args()

    # read tab-seperated parameters file
    parameters = get_model_param_dict(args)
    print(parameters)
    np.random.seed(args.seed)

    # get model architecture to load
    #model, architecture_module=get_model(args, parameters)

    # initialize generators to load data
    train_generator = initializers.initialize_generators(args, "train", parameters, return_coords=True)
    valid_generator = initializers.initialize_generators(args, "valid", parameters, return_coords=True)

    model = load_model_wrapper("/mnt/lab_data2/anusri/chrombpnet/results/tobias/ATAC_PE/GM12878/ATAC_PE_04.07.2022_tobias_corrected_not_softmax/tobias_model/tobias.h5")

    #counts_sum, counts, coords = predict_on_batch_wrapper(valid_generator)
    #print(sum(np.array(list(coords[:,3]))=='0'))
    #print(len(coords[:,3]))
    #print(coords[coords[:,3]=='0'])
    #print(coords[coords[:,0]=='chr19'])

    #print(counts_sum)
    #print(counts)
    train_generator.on_epoch_end()
    #counts_sum, counts, coords, prof_preds, count_preds= predict_on_batch_wrapper(valid_generator)
    #print(sum(np.array(list(coords[:,3]))=='0'))
    #print(len(coords[:,3]))
    #print(coords[coords[:,3]=='0'])
    #print(coords[coords[:,0]=='chr19'])

    #print(counts_sum)
    #print(counts)
    train_generator.on_epoch_end()
    counts_sum, counts, coords, prof_preds, count_preds = predict_on_batch_wrapper(valid_generator, model)
    #print(sum(np.array(list(coords[:,3]))=='0'))
    #print(len(coords[:,3]))
    #print(coords[coords[:,3]=='0'])
    #print(coords[coords[:,0]=='chr19'])

    for idd in range(counts.shape[0])[0:10]:
        plt.figure()
        print(counts_sum.shape)
        print(counts_sum[idd].shape)
        print(np.log(sum(counts_sum[idd])+1))
        print(counts[idd])
        print(prof_preds.shape)
        #plt.plot(prof_preds[idd],label="predictions", alpha=0.5)
        plt.plot(counts_sum[idd]/sum(counts_sum[idd]),label="label", alpha=0.5)
        plt.legend()

        plt.savefig("temp_images/train_"+str(idd)+".png")


    # train the model using the generators
    #fit_and_evaluate(model, train_generator, valid_generator, args, architecture_module)

    # store arguments and and parameters to checkpoint
    #with open(args.output_prefix+'.args.json', 'w') as fp:
    #    json.dump(args.__dict__, fp,  indent=4)
    #with open(args.output_prefix+'.params.json', 'w') as fp:
    #    json.dump(parameters, fp,  indent=4)


if __name__=="__main__":
    main()
