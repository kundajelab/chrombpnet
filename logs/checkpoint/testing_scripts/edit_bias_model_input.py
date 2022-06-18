from tensorflow import keras  # or import keras for standalone version
from tensorflow.keras.layers import Input
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import get_custom_objects
import tensorflow as tf
from tensorflow.keras.models import Model
from src.training.utils import losses

from tensorflow.keras import models
from tensorflow.keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Add, Concatenate, Lambda, Flatten
import numpy as np ;
from tensorflow.keras.backend import int_shape
from tensorflow.keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Flatten
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model
import tensorflow as tf
import random as rn
import os 

import argparse

parser = argparse.ArgumentParser(description='input bias model path')
parser.add_argument("-b", "--bias_path", type=str, required=True, help="bias_model_path")
parser.add_argument("-il", "--inputlen", type=int, required=True, help="new input length")

args = parser.parse_args()

os.environ['PYTHONHASHSEED'] = '0'

def getModelGivenModelOptionsAndWeightInits(pretrained_model, new_input_len):


    #default params (can be overwritten by providing model_params file as input to the training function)
    conv1_kernel_size=21
    profile_kernel_size=75
    num_tasks=1 # not using multi tasking
    
    filters=128
    n_dil_layers=4
    sequence_len=2114
    out_pred_len=1000
    counts_loss_weight=1

    #read in arguments
    seed=1234
    np.random.seed(seed)    
    tf.random.set_seed(seed)
    rn.seed(seed)


    new_input = Input(shape=(new_input_len, 4), name='sequence')
    # Step 1.2 - Crop to match size of the required output size
    cropsize = int(int(new_input_len/2)-int(2114/2))
    inp = Cropping1D(cropsize, name='crop_inp')(new_input)


    # first convolution without dilation
    x = Conv1D(filters,
                kernel_size=conv1_kernel_size,
                padding='valid', 
                activation='relu',
                name='bpnet_1st_conv')(inp)

    layer_names = [str(i) for i in range(1,n_dil_layers+1)]
    for i in range(1, n_dil_layers + 1):
        # dilated convolution
        conv_layer_name = 'bpnet_{}conv'.format(layer_names[i-1])
        conv_x = Conv1D(filters, 
                        kernel_size=3, 
                        padding='valid',
                        activation='relu', 
                        dilation_rate=2**i,
                        name=conv_layer_name)(x)

        x_len = int_shape(x)[1]
        conv_x_len = int_shape(conv_x)[1]
        assert((x_len - conv_x_len) % 2 == 0) # Necessary for symmetric cropping

        x = Cropping1D((x_len - conv_x_len) // 2, name="bpnet_{}crop".format(layer_names[i-1]))(x)
        x = add([conv_x, x])

    # Branch 1. Profile prediction
    # Step 1.1 - 1D convolution with a very large kernel
    prof_out_precrop = Conv1D(filters=num_tasks,
                        kernel_size=profile_kernel_size,
                        padding='valid',
                        name='prof_out_precrop')(x)

    # Step 1.2 - Crop to match size of the required output size
    cropsize = int(int_shape(prof_out_precrop)[1]/2)-int(out_pred_len/2)
    assert cropsize>=0
    assert (cropsize % 2 == 0) # Necessary for symmetric cropping
    prof = Cropping1D(cropsize,
                name='logits_profile_predictions_preflatten')(prof_out_precrop)

    # Branch 2. Counts prediction
    # Step 2.1 - Global average pooling along the "length", the result
    #            size is same as "filters" parameter to the BPNet function

    profile_out = Flatten(name="logits_profile_predictions")(prof)

    gap_combined_conv = GlobalAvgPool1D(name='gap')(x) # acronym - gapcc

    # Step 2.3 Dense layer to predict final counts
    count_out = Dense(num_tasks, name="logcount_predictions")(gap_combined_conv)


    # instantiate keras Model with inputs and outputs
    model=Model(inputs=[new_input],outputs=[profile_out, count_out])

    for layer in model.layers:
        try: 

            #print("copied",layer.name)
            #print("before", layer.get_weights())
            layer.set_weights(pretrained_model.get_layer(name=layer.name).get_weights())    
            #print("after", layer.get_weights())
 
        except:
            print("did not copy",layer.name)
            pass

    model.compile(optimizer=Adam(learning_rate=0.001),
                    loss=[losses.multinomial_nll,'mse'],
                    loss_weights=[1,counts_loss_weight])

    return model 

# or kerassurgeon for standalone Keras
#from tfkerassurgeon import delete_layer, insert_layer


def change_model(model, new_input_shape=(None, 4114, 4),custom_objects=None):
    # replace input shape of first layer


    new_model.add(ResNet)
    

    config = model.layers[0].get_config()
    print(config)
    config['batch_input_shape']=new_input_shape
    model._layers[0]=model.layers[0].from_config(config)

    # rebuild model architecture by exporting and importing via json
    new_model = tf.keras.models.model_from_json(model.to_json(),custom_objects=custom_objects)

    # copy weights from old model to new one
    for layer in new_model._layers:
        try:
            layer.set_weights(model.get_layer(name=layer.name).get_weights())
            print("Loaded layer {}".format(layer.name))
        except:
            print("Could not transfer weights for layer {}".format(layer.name))

    return new_model


def load_model_wrapper(model_h5):
    # read .h5 model
    custom_objects={"tf": tf, "multinomial_nll":losses.multinomial_nll}
    #custom_objects={"tf": tf}
    get_custom_objects().update(custom_objects)    
    model=load_model(model_h5, compile=False)
    print("got the model")
    #model.summary()
    return model, custom_objects

model, custom_objects = load_model_wrapper(args.bias_path+"/bias.h5")

new_model = getModelGivenModelOptionsAndWeightInits(model, args.inputlen)

print(new_model.summary())

new_model.save(args.bias_path+"/bias_input_"+str(args.inputlen)+".h5")
