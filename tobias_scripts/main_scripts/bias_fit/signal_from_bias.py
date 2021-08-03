import pdb 
import numpy as np ;
import tensorflow as tf
import random as rn
import os

os.environ['PYTHONHASHSEED'] = '0'

from keras.backend import int_shape
from sklearn.metrics import average_precision_score
from kerasAC.metrics import * 
from kerasAC.custom_losses import *

import keras;

#import the various keras layers 
from keras.layers import Dense,Activation,Dropout,Flatten,Reshape,Input, Concatenate, Cropping1D, Add
from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers.convolutional import Conv1D
from keras.layers.pooling import GlobalMaxPooling1D,MaxPooling1D,GlobalAveragePooling1D
from keras.layers.normalization import BatchNormalization

from keras.optimizers import Adam
from keras.constraints import maxnorm;
from keras.regularizers import l1, l2    

from keras.models import Model
from tensorflow.keras.models import load_model, model_from_json

def load_model_weights(weight_file,model):
    if weight_file is None:
        #nothing to do
        return model
    #sync the model locally if it's on AWS
    if weight_file.startswith("s3://"):
        s3_model_weights=download_s3_file(weight_file)
    else:
        s3_model_weights=weight_file
    import  h5py
    try:
        model.load_weights(s3_model_weights)#,by_name=True)
    except:
        with h5py.File(s3_model_weights,'r') as file:
            weight_file=file['model_1']
            for layer in model.layers:
                try:
                    layer_weights=weight_file[layer.name]
                except:
                    print('no weights files saved for layer:'+str(layer.name))
                    continue
                try:
                    weights = []
                    # Extract weights
                    for term in layer_weights:
                        if isinstance(layer_weights[term], h5py.Dataset):
                            # Convert weights to numpy array and prepend to list
                            weights.insert(0, np.array(layer_weights[term]))
                    # Load weights to model
                    layer.set_weights(weights)
                    print("loaded weights for layer:"+str(layer.name))
                except Exception as e:
                    print("Error: Could not load weights for layer:"+str(layer.name))
                    raise e
    return model

def load_pretrained_bias(json_string, weights, model_hdf5=None):
    from keras.models import load_model
    from keras.utils.generic_utils import get_custom_objects
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
    if model_hdf5:
        model=load_model(args.model_hdf5)
        print(model.summary())
    else:
        model=model_from_json(open(json_string,"r").read())
        model=load_model_weights(weights,model)

    pretrained_bias_model=model

    #freeze the model
    num_layers=len(pretrained_bias_model.layers)
    for i in range(num_layers):
        pretrained_bias_model.layers[i].trainable=False 
    return pretrained_bias_model 



def get_model_param_dict(param_file):
    '''
    param_file has 2 columns -- param name in column 1, and param value in column 2
    '''
    params={}
    if param_file is None:
        return  params
    for line in open(param_file,'r').read().strip().split('\n'):
        tokens=line.split('\t')
        params[tokens[0]]=tokens[1]
    return params 

def getModelGivenModelOptionsAndWeightInits(args):
    #default params (can be overwritten by providing model_params file as input to the training function)
    model_params=get_model_param_dict(args.model_params)
    profile_loss_weight=float(model_params['profile_loss_weight'])
    counts_loss_weight=float(model_params['counts_loss_weight'])
    pretrained_bias_model=load_pretrained_bias(model_params['json_string'], model_params["weights"])
    print("loaded pre-trained bias model with frozen layers")
    
    #read in arguments
    seed=int(args.seed)
    init_weights=args.init_weights 
    sequence_flank=int(args.tdb_input_flank[0])
    num_tasks=int(args.num_tasks)

    np.random.seed(seed)
    tf.random.set_seed(seed)
    rn.seed(seed)
    
    seq_len=2*sequence_flank
    out_flank=int(args.tdb_output_flank[0])
    out_pred_len=2*out_flank
    print(seq_len)
    print(out_pred_len)

    #define inputs
    bias_counts_input=Input(shape=(1,),name='control_logcount')
    bias_profile_input=Input(shape=(out_pred_len,1),name='control_profile')

    # conv layer without activation 
    profile_out = Conv1D(1,
                         kernel_size=20,
                         padding='same',
                         activation=None,
                         name='profile_out')(bias_profile_input)
    count_out=Dense(1,activation=None,name='count_out')(bias_counts_input)
    model=Model(inputs=[inp],outputs=[profile_out,
                                     count_out])
    print("got model") 
    model.compile(optimizer=Adam(),
                    loss=[MultichannelMultinomialNLL(1),'mse'],
                    loss_weights=[profile_loss_weight,counts_loss_weight])
    print("compiled model")
    return model 


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser(description="view model arch")
    parser.add_argument("--seed",type=int,default=1234)
    parser.add_argument("--init_weights",default=None)
    parser.add_argument("--tdb_input_flank",nargs="+",default=[673])
    parser.add_argument("--tdb_output_flank",nargs="+",default=[500])
    parser.add_argument("--num_tasks",type=int,default=1)
    parser.add_argument("--model_params",default=None)
    
    args=parser.parse_args()
    model=getModelGivenModelOptionsAndWeightInits(args)
    print(model.summary())
    pdb.set_trace() 
    
