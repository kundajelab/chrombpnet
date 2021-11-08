from __future__ import division, print_function, absolute_import
import importlib
import imp
import tensorflow.keras.callbacks as tfcallbacks 
import utils.argmanager as argmanager
import utils.losses as losses
import utils.callbacks as callbacks
import generators.tiledb_generator as tiledb_generator
import generators.batchgen_generator as batchgen_generator
import utils.batchgen_generator_utils as data_utils
import pandas as pd
import splits

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

def initialize_generators(args, parameters, mode, generator_name):

    assert(args.inputlen%2==0)
    assert(args.outputlen%2==0)

    #defaults
    bed_regions=None
    nonpeak_regions=None
    add_revcomp=True
    cts_sum_min_thresh=parameters["cts_sum_min_thresh"]
    cts_sum_max_thresh=parameters["cts_sum_max_thresh"]
    max_jitter=args.max_jitter
    negative_sampling_ratio=args.negative_sampling_ratio

    # get only those peak/non peak regions corresponding to train/valid/test set
    if args.peaks is not None:
            peak_regions=pd.read_csv(args.peaks,header=None,sep='\t',names=NARROWPEAK_SCHEMA)
            peak_regions, chroms=splits.get_bed_regions_for_fold_split(peak_regions, args.fold, mode)

    if args.nonpeaks is not None:
            nonpeak_regions=pd.read_csv(args.nonpeaks,header=None,sep='\t',names=NARROWPEAK_SCHEMA)
            nonpeak_regions, chroms=splits.get_bed_regions_for_fold_split(nonpeak_regions, args.fold, mode) 
    
    if mode=="valid":
            # make sure there is no jitter at test/validation time time - we are testing only at summits
            # make sure there is no reverse complementation
            # fixed negative test for validation and testing
            nonpeak_regions=nonpeak_regions.sample(frac=args.negative_sampling_ratio, replace=False, random_state=args.seed)
            max_jitter=0
            add_revcomp=False
            cts_sum_min_thresh="None" 
            cts_sum_max_thresh="None"
            negative_sampling_ratio=1.0

  
    if generator_name=="tiledb":
        generator=tiledb_generator.TiledbGenerator(
                                        peak_regions=peak_regions,
                                        nonpeak_regions=nonpeak_regions,
                                        genome_fasta=args.genome,
                                        batch_size=args.batch_size,
                                        input_len=args.inputlen,      
                                        output_len=args.outputlen,
                                        max_jitter=max_jitter,
                                        negative_sampling_ratio=args.negative_sampling_ratio,
                                        cts_sum_min_thresh=cts_sum_min_thresh,
                                        cts_sum_max_thresh=cts_sum_max_thresh,
                                        tdb_array=args.tdb_array,
                                        tdb_dataset=args.tdb_dataset,
                                        seed=args.seed,
                                        add_revcomp=add_revcomp,
                                        return_coords=False
                                        )

    if generator_name=="batchgen":
        generator=batchgen_generator.ChromBPNetBatchGenerator(
                                        peak_regions=peak_regions,
                                        nonpeak_regions=nonpeak_regions,
                                        genome_fasta=args.genome,
                                        batch_size=args.batch_size,
                                        inputlen=args.inputlen,                                        
                                        outputlen=args.outputlen,
                                        max_jitter=max_jitter,
                                        negative_sampling_ratio=args.negative_sampling_ratio,
                                        cts_sum_min_thresh=cts_sum_min_thresh,
                                        cts_sum_max_thresh=cts_sum_max_thresh,
                                        cts_bw_file=args.bigwig,
                                        seed=args.seed,
                                        add_revcomp=add_revcomp,
                                        return_coords=False
                                        )
    
    print("data generator is ready!")

    return generator
  

def get_model(args):
    # read model from file - can read model from pre-trained too - but not needed now - will add later
    custom_objects={"MultichannelMultinomialNLL": losses.MultichannelMultinomialNLL}        
    architecture_module=imp.load_source('',args.architecture_from_file)
    model=architecture_module.getModelGivenModelOptionsAndWeightInits(args)
    print("got the model")
    model.summary()
    return model


def fit_and_evaluate(model,train_gen,valid_gen,args):
    model_output_path_string = args.output_prefix
    model_output_path_tf_name=model_output_path_string+".tf"
    model_output_path_logs_name=model_output_path_string+".log"
    model_output_path_arch_name=model_output_path_string+".arch"
    model_output_path_weights_name=model_output_path_string+".weights"

    checkpointer = tfcallbacks.ModelCheckpoint(filepath=model_output_path_tf_name, verbose=1, save_best_only=True, save_weights_only=False)
    earlystopper = tfcallbacks.EarlyStopping(monitor='val_loss', patience=args.early_stop, verbose=1, restore_best_weights=True)
    history= callbacks.LossHistory(model_output_path_logs_name+".batch",args.trackables)
    csvlogger = tfcallbacks.CSVLogger(model_output_path_logs_name, append = False)
    reduce_lr = tfcallbacks.ReduceLROnPlateau(monitor='val_loss', factor=0.4,patience=args.early_stop-2, min_lr=0.00000001)
    cur_callbacks=[checkpointer,earlystopper,csvlogger,reduce_lr,history]

    model.fit(train_gen,
              validation_data=valid_gen,
              epochs=args.epochs,
              verbose=1,
              callbacks=cur_callbacks,
              shuffle=True)

    # save model as hdf5 too - add this - ##
    print('fit_generator complete') 
    model.save_weights(model_output_path_weights_name)

    print('weights saved') 
    architecture_string=model.to_json()
    with open(model_output_path_arch_name,'w') as outf:
        outf.write(architecture_string)
    print('saved model architecture')


def main():

    # read arguments
    args=argmanager.fetch_train_chrombpnet_args()

    # read tab-seperated paprmeters file
    data = open(args.params).readlines()
    parameters={}
    for line in data:
        [key,value] = line.strip().split("\t")
        parameters[key] = value

    print(parameters)
    # initialize generators to load data
    train_generator = initialize_generators(args, parameters, mode="train", generator_name=args.generator)
    valid_generator = initialize_generators(args, parameters, mode="valid", generator_name=args.generator)

    # get model architecture to load
    model=get_model(args)

    # train the model using the generators
    fit_and_evaluate(model, train_generator, valid_generator, args)
    

if __name__=="__main__":
    main()
