from __future__ import division, print_function, absolute_import
import importlib.machinery
import tensorflow.keras.callbacks as tfcallbacks 
import chrombpnet.training.utils.argmanager as argmanager
import chrombpnet.training.utils.losses as losses
import chrombpnet.training.utils.callbacks as callbacks
import chrombpnet.training.data_generators.initializers as initializers
import pandas as pd
import os
import json
import numpy as np

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
os.environ['PYTHONHASHSEED'] = '0'

def get_model(args, parameters):
    """
    Read a model definition from a python file. This function can be used to read any model architecture that takes sequence as input
    and outputs a two task model.  Task one to predict the probability distribution of a profile and task two to predict the total counts in a profile.
    Look at .py models in src/training/models/ for examples. I will try to provide a dummy model as example here - for later.
    The files should have the following two functions - getModelGivenModelOptionsAndWeightInits and save_model_without_bias
    """
    architecture_module=importlib.machinery.SourceFileLoader('',args.architecture_from_file).load_module()
    model=architecture_module.getModelGivenModelOptionsAndWeightInits(args, parameters)
    print("got the model")
    return model, architecture_module

def fit_and_evaluate(model,train_gen,valid_gen,args,architecture_module):
    model_output_path_h5_name=args.output_prefix+".h5"
    model_output_path_logs_name=args.output_prefix+".log"

    checkpointer = tfcallbacks.ModelCheckpoint(filepath=model_output_path_h5_name, monitor="val_loss", mode="min",  verbose=1, save_best_only=True)
    earlystopper = tfcallbacks.EarlyStopping(monitor='val_loss', mode="min", patience=args.early_stop, verbose=1, restore_best_weights=True)
    history= callbacks.LossHistory(model_output_path_logs_name+".batch",args.trackables)
    csvlogger = tfcallbacks.CSVLogger(model_output_path_logs_name, append=False)
    #reduce_lr = tfcallbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.4, patience=args.early_stop-2, min_lr=0.00000001)
    cur_callbacks=[checkpointer,earlystopper,csvlogger,history]

    model.fit(train_gen,
              validation_data=valid_gen,
              epochs=args.epochs,
              verbose=1,
              callbacks=cur_callbacks)

    print('save model') 
    model.save(model_output_path_h5_name)

    architecture_module.save_model_without_bias(model, args.output_prefix)


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
    assert(args.chr_fold_path==params["chr_fold_path"]) # the parameters were generated on a different folds compared to the given fold

    assert(int(params["inputlen"])%2==0)
    assert(int(params["outputlen"])%2==0)

    return params 

def main(args):


    # read tab-seperated parameters file
    parameters = get_model_param_dict(args)
    print(parameters)
    np.random.seed(args.seed)

    # get model architecture to load
    model, architecture_module=get_model(args, parameters)

    # initialize generators to load data
    train_generator = initializers.initialize_generators(args, "train", parameters, return_coords=False)
    valid_generator = initializers.initialize_generators(args, "valid", parameters, return_coords=False)

    # train the model using the generators
    fit_and_evaluate(model, train_generator, valid_generator, args, architecture_module)

    # store arguments and and parameters to checkpoint
    with open(args.output_prefix+'.args.json', 'w') as fp:
        json.dump(args.__dict__, fp,  indent=4)
    #with open(args.output_prefix+'.params.json', 'w') as fp:
    #    json.dump(parameters, fp,  indent=4)


if __name__=="__main__":
    # read arguments
    args=argmanager.fetch_train_args()
    main(args)

