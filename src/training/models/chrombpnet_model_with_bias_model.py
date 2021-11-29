import numpy as np ;
from keras.backend import int_shape
from keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Add, Concatenate, Lambda
from keras.optimizers import Adam
from keras.models import Model
from utils.losses import MultichannelMultinomialNLL
import tensorflow as tf
import random as rn
import os 

os.environ['PYTHONHASHSEED'] = '0'


def load_pretrained_bias(model_hdf5):
    from keras.models import load_model
    from keras.utils.generic_utils import get_custom_objects
    custom_objects={"MultichannelMultinomialNLL": MultichannelMultinomialNLL}    
    get_custom_objects().update(custom_objects)
    pretrained_bias_model=load_model(model_hdf5)
    #freeze the model
    num_layers=len(pretrained_bias_model.layers)
    for i in range(num_layers):
        pretrained_bias_model.layers[i].trainable=False
    return pretrained_bias_model


def getModelGivenModelOptionsAndWeightInits(args, model_params):
    #default params (can be overwritten by providing model_params file as input to the training function)
    filters=512
    n_dil_layers=8
    conv1_kernel_size=21
    profile_kernel_size=75
    counts_loss_weight=1
    num_tasks=1 # not using multi tasking
    
    if 'filters' in model_params:
        filters=int(model_params['filters'])
    if 'n_dil_layers' in model_params:
        n_dil_layers=int(model_params['n_dil_layers'])
    if 'conv1_kernel_size' in model_params:
        conv1_kernel_size=int(model_params['conv1_kernel_size'])
    if 'profile_kernel_size' in model_params:
        profile_kernel_size=int(model_params['profile_kernel_size'])
    if 'cnts_loss_weight' in model_params:
        counts_loss_weight=float(model_params['counts_loss_weight'])
    if 'bias_model_path' in model_params:
        bias_model_path=model_params['bias_model_path']
    else:
        print("Bias model hdf5 path expected")

    bias_model = load_pretrained_bias(bias_model_path)

    print("params:")
    print("filters:"+str(filters))
    print("n_dil_layers:"+str(n_dil_layers))
    print("conv1_kernel_size:"+str(conv1_kernel_size))
    print("profile_kernel_size:"+str(profile_kernel_size))
    print("counts_loss_weight:"+str(counts_loss_weight))
    
    #read in arguments
    seed=args.seed
    np.random.seed(seed)    
    tf.random.set_seed(seed)
    rn.seed(seed)
    sequence_len=int(args.inputlen)
    out_pred_len=int(args.outputlen)

    #define inputs
    inp = Input(shape=(sequence_len, 4),name='sequence')    

    # first convolution without dilation
    x = Conv1D(filters,
                kernel_size=conv1_kernel_size,
                padding='valid', 
                activation='relu',
                name='chrombpnet_1st_conv')(inp)

    layer_names = [str(i) for i in range(1,n_dil_layers+1)]
    for i in range(1, n_dil_layers + 1):
        # dilated convolution
        conv_layer_name = 'chrombpnet_{}conv'.format(layer_names[i-1])
        conv_x = Conv1D(filters, 
                        kernel_size=3, 
                        padding='valid',
                        activation='relu', 
                        dilation_rate=2**i,
                        name=conv_layer_name)(x)

        x_len = int_shape(x)[1]
        conv_x_len = int_shape(conv_x)[1]
        assert((x_len - conv_x_len) % 2 == 0) # Necessary for symmetric cropping

        x = Cropping1D((x_len - conv_x_len) // 2, name="chrombpnet_{}crop".format(layer_names[i-1]))(x)
        x = add([conv_x, x])

    # Branch 1. Profile prediction
    # Step 1.1 - 1D convolution with a very large kernel
    prof_out_precrop = Conv1D(filters=num_tasks,
                        kernel_size=profile_kernel_size,
                        padding='valid',
                        name='chrombpnet_prof_out_precrop')(x)

    # Step 1.2 - Crop to match size of the required output size
    cropsize = int(int_shape(prof_out_precrop)[1]/2)-int(out_pred_len/2)
    assert cropsize>=0
    assert (cropsize % 2 == 0) # Necessary for symmetric cropping
    profile_out_prebias = Cropping1D(cropsize,
                name='chrombpnet_logits_profile_predictions')(prof_out_precrop)

    # Branch 2. Counts prediction
    # Step 2.1 - Global average pooling along the "length", the result
    #            size is same as "filters" parameter to the BPNet function
    gap_combined_conv = GlobalAvgPool1D(name='gap')(x) # acronym - gapcc

    # Step 2.3 Dense layer to predict final counts
    count_out_prebias = Dense(num_tasks, name="chrombpnet_logcount_predictions")(gap_combined_conv)

    ## get bias output
    bias_output=bias_model(inp)
    profile_out = Add(name="logits_profile_predictions")([profile_out_prebias,bias_output[0]])
    concat_counts = Concatenate(axis=-1)([count_out_prebias, bias_output[1]])
    count_out = Lambda(lambda x: tf.math.reduce_logsumexp(x, axis=-1, keepdims=True),
                        name="logcount_predictions")(concat_counts)

    # instantiate keras Model with inputs and outputs
    model=Model(inputs=[inp],outputs=[profile_out, count_out])

    model.compile(optimizer=Adam(learning_rate=args.learning_rate),
                    loss=[MultichannelMultinomialNLL(num_tasks),'mse'],
                    loss_weights=[1,counts_loss_weight])

    return model 


def get_model_without_bias(model, output_prefix):
    profile_output_without_bias = model.get_layer("chrombpnet_logits_profile_predictions").output
    counts_output_without_bias = model.get_layer("chrombpnet_logcount_predictions").output
    model_without_bias = Model(inputs=model.inputs,outputs=[profile_output_without_bias, counts_output_without_bias])
    print('save model without bias') 
    model_without_bias.save(args.output_prefix+"wo_bias.h5")