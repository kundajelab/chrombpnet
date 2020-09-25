import pdb 
import numpy as np ;
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

def load_pretrained_bias(model_hdf5):
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
    return pretrained_bias_model 


def getModelGivenModelOptionsAndWeightInits(args):
    #default params (can be overwritten by providing model_params file as input to the training function)
    filters=500
    n_dil_layers=6
    conv1_kernel_size=21
    profile_kernel_size=75
    control_smoothing=[1, 50]
    counts_loss_weight=1
    profile_loss_weight=1
    
    model_params=get_model_param_dict(args.model_params)
    if 'filters' in model_params:
        filters=int(model_params['filters'])
    if 'n_dil_layers' in model_params:
        n_dil_layers=int(model_params['n_dil_layers'])
    if 'conv1_kernel_size' in model_params:
        conv1_kernel_size=int(model_params['conv1_kernel_size'])
    if 'profile_kernel_size' in model_params:
        profile_kernel_size=int(model_params['profile_kernel_size'])
    if 'control_smoothing' in model_params:
        control_smoothing=[int(i) for i in model_params['control_smoothing'].split(',')]
    if 'counts_loss_weight' in model_params:
        counts_loss_weight=float(model_params['counts_loss_weight'])
    if 'profile_loss_weight' in model_params:
        profile_loss_weight=float(model_params['profile_loss_weight'])

    print("params:")
    print("filters:"+str(filters))
    print("n_dil_layers:"+str(n_dil_layers))
    print("conv1_kernel_size:"+str(conv1_kernel_size))
    print("profile_kernel_size:"+str(profile_kernel_size))
    print("control_smoothing:"+str(control_smoothing))
    print("counts_loss_weight:"+str(counts_loss_weight))
    print("profile_loss_weight:"+str(profile_loss_weight))
    
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

    #load the pretrained bias model
    pretrained_bias_model=load_pretrained_bias(model_params['pretrained_bias_model'])
    
    #define inputs
    inp = Input(shape=(seq_len, 4),name='sequence')    

    # first convolution without dilation
    first_conv = Conv1D(filters,
                        kernel_size=conv1_kernel_size,
                        padding='valid', 
                        activation='relu',
                        name='1st_conv')(inp)
    # 6 dilated convolutions with resnet-style additions
    # each layer receives the sum of feature maps 
    # from all previous layers
    res_layers = [(first_conv, '0_non_dil')] 
    layer_names = [str(i)+"_dil" for i in range(1,n_dil_layers+1)]
    for i in range(1, n_dil_layers + 1):
        if i == 1:
            res_layers_sum = first_conv
        else:
            res_layers_sum = Add(name='add'+str(i))([l for l, _ in res_layers])

        # dilated convolution
        conv_layer_name = '{}conv'.format(layer_names[i-1])
        conv_output = Conv1D(filters, 
                             kernel_size=3, 
                             padding='valid',
                             activation='relu', 
                             dilation_rate=2**i,
                             name=conv_layer_name)(res_layers_sum)

        # get shape of latest layer and crop 
        # all other previous layers in the list to that size
        conv_output_shape =int_shape(conv_output)
        cropped_layers = []
        for lyr, name in res_layers:
            lyr_shape =int_shape(lyr)
            cropsize = int(lyr_shape[1]/2) - int(conv_output_shape[1]/2)
            lyr_name = '{}-crop_{}th_dconv'.format(name.split('-')[0], i)
            cropped_layers.append((Cropping1D(cropsize,
                                              name=lyr_name)(lyr),
                                  lyr_name))
        
        # append to the list of previous layers
        cropped_layers.append((conv_output, conv_layer_name))
        res_layers = cropped_layers

    # the final output from the 6 dilated convolutions 
    # with resnet-style connections
    combined_conv = Add(name='combined_conv')([l for l, _ in res_layers])
    profile_out_prebias = Conv1D(filters=num_tasks,
                                 kernel_size=profile_kernel_size,
                                 padding='valid',
                                 name='profile_out_prebias')(combined_conv)
    # Step 1.2 - Crop to match size of the required output size, a minimum
    #            difference of 346 is required between input seq len and ouput len
    profile_out_prebias_shape =int_shape(profile_out_prebias)
    cropsize = int(profile_out_prebias_shape[1]/2)-int(out_pred_len/2)
    profile_out_prebias = Cropping1D(cropsize,
                                     name='prof_out_crop2match_output')(profile_out_prebias)

    #ADD IN THE BIAS
    bias_output=pretrained_bias_model(inp)
    #concatenate bias output with profile_out_prebias
    concat_with_bias_prof=Concatenate(axis=-1)([profile_out_prebias,bias_output[0]])
    profile_out = Conv1D(filters=num_tasks,
                         kernel_size=1,
                         name="profile_predictions")(concat_with_bias_prof)
    #COUNTS ARM 
    gap_combined_conv = GlobalAveragePooling1D(name='gap')(combined_conv)
    #concatenate gap_combined_conv with bias output
    concat_with_bias_count=Concatenate(axis=-1)([gap_combined_conv,bias_output[1]])
    count_out = Dense(num_tasks, name="logcount_predictions")(concat_with_bias_count)
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
    
