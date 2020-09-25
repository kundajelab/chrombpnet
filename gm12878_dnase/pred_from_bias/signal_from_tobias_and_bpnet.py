import pdb 
import numpy as np
from kerasAC.metrics import * 
from kerasAC.custom_losses import *
import keras

#import the various keras layers 
from keras.layers import Dense,Activation,Input, Concatenate, Cropping1D, Add
from keras.layers.core import Dropout, Reshape, Dense, Activation, Flatten
from keras.layers.convolutional import Conv1D, Conv2D
from keras.layers.pooling import GlobalMaxPooling1D,MaxPooling1D,GlobalAveragePooling1D
from keras.layers.normalization import BatchNormalization

from keras.optimizers import Adam
from keras.constraints import maxnorm
from keras.regularizers import l1, l2    
from keras.models import Model


def load_pretrained_bias(model_hdf5,model_name=None):
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
    pretrained_bias_model=load_model(model_hdf5)
    #freeze the model
    num_layers=len(pretrained_bias_model.layers)
    for i in range(num_layers):
        pretrained_bias_model.layers[i].trainable=False
    if model_name is not None:
        pretrained_bias_model.name=model_name 
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
    pretrained_model_names=model_params['pretrained_bias_model'].split(',')
    assert len(pretrained_model_names)==2
    pretrained_bias_models=[load_pretrained_bias(pretrained_model_names[i],'model_'+str(i)) for i in range(len(pretrained_model_names))]
    print("loaded pre-trained bias model with frozen layers")
    
    #read in arguments
    seed=args.seed
    init_weights=args.init_weights 
    sequence_flank=args.tdb_input_flank[0]
    num_tasks=args.num_tasks
    
    seq_len=2*sequence_flank
    out_flank=args.tdb_output_flank[0]
    out_pred_len=2*out_flank
    print(seq_len)
    print(out_pred_len)
    #define inputs
    inp = Input(shape=(seq_len, 4),name='sequence')    
    bias_output0 = pretrained_bias_models[0] (inp)
    bias_output1 = pretrained_bias_models[1] (inp)
    #stack the outputs
    concat_biases_profile=Concatenate(axis=-1,name='concat_biases_profile')([bias_output0[0],bias_output1[0]])
    concat_biases_count=Concatenate(axis=-1,name='concat_biases_count')([bias_output0[1],bias_output1[1]])

    print("concat_biases_profile shape:"+str(concat_biases_profile))
    print("concat_biases_count shape:"+str(concat_biases_count))
    
    # conv layer without activation 
    profile_out = Conv1D(1,
                         kernel_size=20,
                         padding='same',
                         activation=None,
                         name='profile_out')(concat_biases_profile)
    
    count_out=Dense(1,activation=None,name='count_out')(concat_biases_count)
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
    
