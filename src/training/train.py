from __future__ import division, print_function, absolute_import
import importlib.machinery
import tensorflow.keras.callbacks as tfcallbacks 
import utils.argmanager as argmanager
import utils.losses as losses
import utils.callbacks as callbacks
import generators.initializers
import pandas as pd
import splits
import os
import json

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
os.environ['PYTHONHASHSEED'] = '0'

def get_model(args, parameters):
    # read model from file - can read model from pre-trained too - but not needed now - will add later
    architecture_module=importlib.machinery.SourceFileLoader('',args.architecture_from_file).load_module()
    model=architecture_module.getModelGivenModelOptionsAndWeightInits(args, parameters)
    print("got the model")
    model.summary()
    return model

def fit_and_evaluate(model,train_gen,valid_gen,args):
    model_output_path_string = args.output_prefix
    model_output_path_h5_name=model_output_path_string+".h5"
    model_output_path_logs_name=model_output_path_string+".log"
    model_output_path_arch_name=model_output_path_string+".arch"
    model_output_path_weights_name=model_output_path_string+".weights"

    checkpointer = tfcallbacks.ModelCheckpoint(filepath=model_output_path_h5_name, monitor="val_loss", mode="min",  verbose=1, save_best_only=True)
    earlystopper = tfcallbacks.EarlyStopping(monitor='val_loss', mode="min", patience=args.early_stop, verbose=1, restore_best_weights=True)
    history= callbacks.LossHistory(model_output_path_logs_name+".batch",args.trackables)
    csvlogger = tfcallbacks.CSVLogger(model_output_path_logs_name, append=False)
    reduce_lr = tfcallbacks.ReduceLROnPlateau(monitor='val_loss',factor=0.4, patience=args.early_stop-2, min_lr=0.00000001)
    cur_callbacks=[checkpointer,earlystopper,csvlogger,reduce_lr,history]

    model.fit(train_gen,
              validation_data=valid_gen,
              epochs=args.epochs,
              verbose=1,
              max_queue_size=100,
              workers=10,
              callbacks=cur_callbacks,
              shuffle=True)

    print('fit_generator complete') 
    print('save model') 
    model.save(model_output_path_h5_name)

    model.save_weights(model_output_path_weights_name)

    architecture_string=model.to_json()
    with open(model_output_path_arch_name,'w') as outf:
        outf.write(architecture_string)

def get_model_param_dict(param_file):
    '''
    param_file has 2 columns -- param name in column 1, and param value in column 2
    You can pass model specfic parameters to design your own model with this
    '''
    params={}
    if param_file is None:
        return  params
    for line in open(param_file,'r').read().strip().split('\n'):
        tokens=line.split('\t')
        params[tokens[0]]=tokens[1]
    return params 

def main():

    # read arguments
    args=argmanager.fetch_train_chrombpnet_args()

    # read tab-seperated paprmeters file
    parameters = get_model_param_dict(args.params)
    print(parameters)

    # get model architecture to load
    model=get_model(args, parameters)

    # initialize generators to load data
    train_generator = initializers.initialize_generators(args, "train", args.generator, parameters, False)
    valid_generator = initializers.initialize_generators(args, "valid", args.generator, parameters, False)

    # train the model using the generators
    fit_and_evaluate(model, train_generator, valid_generator, args)

    # store args and and store params to checkpoint
    with open(args.output_prefix+'.args.json', 'w') as fp:
        json.dump(args.__dict__, fp,  indent=4)
    with open(args.output_prefix+'.params.json', 'w') as fp:
        json.dump(parameters, fp,  indent=4)


if __name__=="__main__":
    main()
