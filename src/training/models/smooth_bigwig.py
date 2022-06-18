import numpy as np ;
from tensorflow.keras.backend import int_shape
from tensorflow.keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Flatten
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.models import Model
from utils.losses import multinomial_nll
import tensorflow as tf
import random as rn
import os 

os.environ['PYTHONHASHSEED'] = '0'

def getModelGivenModelOptionsAndWeightInits(args, model_params):
    #default params (can be overwritten by providing model_params file as input to the training function)
    num_tasks=1 # not using multi tasking
    
    filters=int(model_params['filters'])
    n_dil_layers=int(model_params['n_dil_layers'])
    counts_loss_weight=float(model_params['counts_loss_weight'])
    sequence_len=int(model_params["inputlen"])
    out_pred_len=int(model_params["outputlen"])

    print("params:")
    print("filters:"+str(filters))
    print("n_dil_layers:"+str(n_dil_layers))
    print("counts_loss_weight:"+str(counts_loss_weight))
    
    #read in arguments
    seed=args.seed
    np.random.seed(seed)    
    tf.random.set_seed(seed)
    rn.seed(seed)

    #define inputs
    inp_bias_logits = Input(shape=(out_pred_len+100,1), name="bias_logits")
    inp_bias_logcounts = Input(shape=(1,), name="bias_logcounts")

    # first convolution without dilation
    x = Conv1D(1,
                kernel_size=75,
                padding='valid', 
                activation=None,
                name='bpnet_1st_conv')(inp_bias_logits)

    # Step 1.2 - Crop to match size of the required output size
    cropsize = int(int_shape(x)[1]/2)-int(out_pred_len/2)
    assert cropsize>=0
    assert (int_shape(x)[1] % 2 == 0)
    #assert (cropsize % 2 == 0) # Necessary for symmetric cropping
    prof = Cropping1D(cropsize,
                name='logits_profile_predictions_preflatten')(x)

    # Branch 2. Counts prediction
    # Step 2.1 - Global average pooling along the "length", the result
    #            size is same as "filters" parameter to the BPNet function

    profile_out = Flatten(name="logits_profile_predictions")(prof)


    # Step 2.3 Dense layer to predict final counts
    count_out = Dense(num_tasks, name="logcount_predictions")(inp_bias_logcounts)

    # instantiate keras Model with inputs and outputs
    model=Model(inputs=[inp_bias_logits, inp_bias_logcounts],outputs=[profile_out, count_out])

    model.compile(optimizer=Adam(learning_rate=args.learning_rate),
                    loss=[multinomial_nll,'mse'],
                    loss_weights=[1,counts_loss_weight])

    return model 

def save_model_without_bias(model, output_prefix):
    # nothing to do 
    # all model architectures have this function
    # defining this tosafeguard if the users uses the arugument save_model_without_bias argument on bias model accidentally 
    return
