import keras
from keras.layers import Input, Cropping1D, add, Conv1D, GlobalAvgPool1D, Dense, Flatten
from keras.models import Model
from chrombpnet.training.utils.losses import multinomial_nll
from chrombpnet.training.optimizers import make_optimizer
from chrombpnet.training.runtime import head_dtype
import os 

os.environ['PYTHONHASHSEED'] = '0'

def getModelGivenModelOptionsAndWeightInits(args, model_params):
    #default params (can be overwritten by providing model_params file as input to the training function)
    conv1_kernel_size=21
    profile_kernel_size=75
    num_tasks=1 # not using multi tasking
    
    filters=int(model_params['filters'])
    n_dil_layers=int(model_params['n_dil_layers'])
    counts_loss_weight=float(model_params['counts_loss_weight'])
    sequence_len=int(model_params["inputlen"])
    out_pred_len=int(model_params["outputlen"])

    print("params:")
    print("filters:"+str(filters))
    print("n_dil_layers:"+str(n_dil_layers))
    print("conv1_kernel_size:"+str(conv1_kernel_size))
    print("profile_kernel_size:"+str(profile_kernel_size))
    print("counts_loss_weight:"+str(counts_loss_weight))
    
    #read in arguments
    seed=args.seed
    keras.utils.set_random_seed(seed)
    # float32 output heads under --precision bf16 (None = the global dtype policy)
    out_dtype=head_dtype()

    #define inputs
    inp = Input(shape=(sequence_len, 4),name='sequence')    

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

        x_len = x.shape[1]
        conv_x_len = conv_x.shape[1]
        assert((x_len - conv_x_len) % 2 == 0) # Necessary for symmetric cropping

        x = Cropping1D((x_len - conv_x_len) // 2, name="bpnet_{}crop".format(layer_names[i-1]))(x)
        x = add([conv_x, x])

    # Branch 1. Profile prediction
    # Step 1.1 - 1D convolution with a very large kernel
    prof_out_precrop = Conv1D(filters=num_tasks,
                        kernel_size=profile_kernel_size,
                        padding='valid',
                        name='prof_out_precrop',
                        dtype=out_dtype)(x)

    # Step 1.2 - Crop to match size of the required output size
    cropsize = int(prof_out_precrop.shape[1]/2)-int(out_pred_len/2)
    assert cropsize>=0
    assert (prof_out_precrop.shape[1] % 2 == 0) # Necessary for symmetric cropping
    prof = Cropping1D(cropsize,
                name='logits_profile_predictions_preflatten',
                dtype=out_dtype)(prof_out_precrop)

    # Branch 2. Counts prediction
    # Step 2.1 - Global average pooling along the "length", the result
    #            size is same as "filters" parameter to the BPNet function

    # the output layer names are a contract: they make the log keys (logits_profile_predictions_loss,
    # logcount_predictions_loss) and find_chrombpnet_hyperparams looks up the count Dense by name
    profile_out = Flatten(name="logits_profile_predictions", dtype=out_dtype)(prof)

    gap_combined_conv = GlobalAvgPool1D(name='gap', dtype=out_dtype)(x) # acronym - gapcc

    # Step 2.3 Dense layer to predict final counts
    count_out = Dense(num_tasks, name="logcount_predictions", dtype=out_dtype)(gap_combined_conv)

    # instantiate keras Model with inputs and outputs
    model=Model(inputs=inp,outputs=[profile_out, count_out])

    model.compile(optimizer=make_optimizer(args),
                    loss=[multinomial_nll,'mse'],
                    loss_weights=[1,counts_loss_weight])

    return model 

def save_model_without_bias(model, output_prefix):
    # nothing to do 
    # all model architectures have this function
    # defining this tosafeguard if the users uses the arugument save_model_without_bias argument on bias model accidentally 
    return
