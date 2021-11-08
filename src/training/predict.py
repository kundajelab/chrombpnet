
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import time
#graceful shutdown
import psutil
import signal 
import os
#multithreading
#from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Pool,Process, Queue 
import warnings
import numpy as np
import pysam
import pandas as pd
import tensorflow as tf 
from kerasAC.activations import softMaxAxis1
from .calibrate import * 
from .generators.basic_generator import *
from .generators.tiledb_predict_generator import *
from .tiledb_config import *
from .s3_sync import *
from .get_model import *
from .splits import *
from kerasAC.config import args_object_from_args_dict
from kerasAC.performance_metrics import *
from kerasAC.custom_losses import *
from kerasAC.metrics import recall, specificity, fpr, fnr, precision, f1
import argparse
import yaml 
import h5py 
import pickle
import numpy as np 
import tensorflow.keras as keras 
from keras.losses import *
from keras.models import Model
#from keras.utils import multi_gpu_model
from kerasAC.custom_losses import *
from abstention.calibration import PlattScaling, IsotonicRegression 
import random
import pdb 

def write_predictions(args):
    '''
    keys lab_i, pred_i indicate the labels and predictions matrices for loss i. 
    The tasks are encoded in the final axis. i.e. profile labels matrix would be (N,L,T), where N = number of samples, L= length of profile, T = number of tasks.
    '''
    try:
        if args.predictions_and_labels_hdf5.startswith('s3://'):
            #use a local version of the file and upload to s3 when finished
            bucket,filename=s3_string_parse(args.predictions_and_labels_hdf5)
            out_predictions_prefix=filename.split('/')[-1]+".predictions"
        else: 
            out_predictions_prefix=args.predictions_and_labels_hdf5+".predictions"
        #create output h5py file
        outf=h5py.File(out_predictions_prefix,'w')
        first=True
        while True:
            cur_out=queue.get()
            if type(cur_out)==str:
                if cur_out=="FINISHED":
                    outf.close()
                    return
            coords=cur_out[0]
            labels=cur_out[1]
            pred=cur_out[2]
            num_datasets=len(pred)
            if first is True:
                #create the predictions & labels datasets
                for i in range(num_datasets):
                    cur_pred_shape=pred[i].shape
                    if len(cur_pred_shape)==1:
                        maxshape=(None,)
                    else:
                        maxshape=(None,)+cur_pred_shape[1::]
                        
                    outf.create_dataset('pred_'+str(i),data=pred[i],maxshape=maxshape)
                    outf.create_dataset('lab_'+str(i),data=labels[i],maxshape=maxshape)
                #create the coord dataset
                coord_shape=coords.shape
                if len(coord_shape)==1:
                    maxshape=(None,)
                else:
                    maxshape=(None,)+coord_shape[1::]
                outf.create_dataset('coords',data=coords,maxshape=maxshape)
                

                first=False
            else:
                #append to hdf5
                
                #append coords
                curshape=len(outf['coords'])
                batchshape=coords.shape[0]
                newshape=curshape+batchshape
                outf['coords'].resize(newshape,axis=0)
                outf['coords'][-1*batchshape::]=coords
                
                #append predictions
                for i in range(num_datasets): 
                    curshape=len(outf['pred_'+str(i)])
                    batchshape=pred[i].shape[0]
                    newshape=curshape+batchshape
                    outf['pred_'+str(i)].resize(newshape,axis=0)
                    outf['pred_'+str(i)][-1*batchshape::]=pred[i]                
                    outf['lab_'+str(i)].resize(newshape,axis=0)
                    outf['lab_'+str(i)][-1*batchshape::]=labels[i]                
    except KeyboardInterrupt:
        #shutdown the pool
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise 
    except Exception as e:
        print(e)
        #shutdown the pool
        # Kill remaining child processes
        kill_child_processes(os.getpid())
        raise e


def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def kill_child_processes(parent_pid, sig=signal.SIGTERM):
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        return
    children = parent.children(recursive=True)
    for process in children:
        process.send_signal(sig)
        
def get_batch_wrapper(idx):
    X,y,coords=test_generator[idx]
    if type(y) is not list:
        y=[y]
    try:
        y=[i.squeeze(axis=-1) for i in y]
    except:
        pass
    if type(X) is not list:
        X=[X]
    return X,y,coords
    #represent coords w/ string, MultiIndex not supported in table append mode
    #set the column names for the MultiIndex
    #coords=pd.MultiIndex.from_tuples(coords,names=['CHR','CENTER','STRAND'])
    #y=[pd.DataFrame(i,index=coords) for i in y]
    #return X,y,coords


def get_tiledb_predict_generator(args):
    global test_generator
    if args.upsample_ratio_list_predict is not None:
        upsample_ratio_predict=args.upsample_ratio_list_predict[0]
        print("warning! only a single ratio for upsampling supported for tiledb as of now")
    else:
        upsample_ratio_predict=None
    import tiledb
    tdb_config=get_default_config() 
    tdb_ctx=tiledb.Ctx(config=tdb_config)
    #you can only specify one (or neither) or args.fold or args.predict chroms 
    assert (args.fold is None) or (args.predict_chroms is None)
    if args.fold is None:
        predict_chroms=args.predict_chroms
    else:
        predict_chroms=get_chroms(args,split='test')
    print(predict_chroms)
#<<<<<<< Updated upstream
    if args.bed_regions is not None:
        bed_regions=pd.read_csv(args.bed_regions,header=None,sep='\t')
        bed_regions_test=get_bed_regions_for_fold_split(bed_regions,args.genome,args.fold,'test')
    else:
        bed_regions_test=args.bed_regions
#=======
    bed_regions_test = args.bed_regions
#>>>>>>> Stashed changes
    test_generator=TiledbPredictGenerator(ref_fasta=args.ref_fasta,
                                          batch_size=args.batch_size,
                                          tdb_array=args.tdb_array,
                                          tdb_partition_attribute_for_upsample=args.tdb_partition_attribute_for_upsample,
                                          tdb_partition_thresh_for_upsample=args.tdb_partition_thresh_for_upsample,
                                          tdb_partition_datasets_for_upsample=args.tdb_partition_datasets_for_upsample,
                                          upsample_ratio=upsample_ratio_predict,
                                          num_threads=args.upsample_threads,
                                          tdb_ambig_attribute=args.tdb_ambig_attribute,
                                          tdb_input_source_attribute=args.tdb_input_source_attribute,
                                          tdb_input_flank=args.tdb_input_flank,
                                          tdb_input_min=args.tdb_input_min,
                                          tdb_input_max=args.tdb_input_max,
                                          tdb_output_source_attribute=args.tdb_output_source_attribute,
                                          tdb_output_flank=args.tdb_output_flank,
                                          tdb_output_min=args.tdb_output_min,
                                          tdb_output_max=args.tdb_output_max,
                                          num_inputs=args.num_inputs,
                                          num_outputs=args.num_outputs,
                                          tdb_input_aggregation=args.tdb_input_aggregation,
                                          tdb_input_transformation=args.tdb_input_transformation,
                                          pseudocount=args.tdb_transformation_pseudocount,
                                          tdb_output_aggregation=args.tdb_output_aggregation,
                                          tdb_output_transformation=args.tdb_output_transformation,                                          
                                          tiledb_stride=args.tiledb_stride,
                                          chrom_sizes=args.chrom_sizes,
                                          chroms=predict_chroms,
                                          tdb_input_datasets=args.tdb_input_datasets,
                                          tdb_output_datasets=args.tdb_output_datasets,
                                          tdb_config=tdb_config,
                                          tdb_ctx=tdb_ctx,
                                          bed_regions=bed_regions_test,
                                          bed_regions_center=args.bed_regions_center,
                                          add_revcomp=args.revcomp)
    print("created TiledbPredictGenerator")    
    return test_generator 

def predict_on_batch_wrapper(args,model,test_generator):
    num_batches=len(test_generator)
    for idx in range(num_batches):
        if idx%100==0:
            print(str(idx)+'/'+str(num_batches))
        X,y,coords=get_batch_wrapper(idx)
        #get the model predictions            
        preds=model.predict_on_batch(X)
        if type(preds) is not list:
            preds=[preds]
        try:
            preds=[i.squeeze(axis=-1) for i in preds]
        except:
            pass 
        queue.put((coords,y,preds))
    print("finished with tiledb predictions!")
    queue.put("FINISHED")
    queue.close() 
    return

def get_model_layer_functor(model,target_layer_idx):
    from keras import backend as K
    inp=model.input
    outputs=model.layers[target_layer_idx].output
    functor=K.function([inp], [outputs])
    return functor 

def get_layer_outputs(functor,X):
    return functor([X])

def predict(args):
    if type(args)==type({}):
        args=args_object_from_args_dict(args) 
    global queue
    
    queue=Queue()
    
    
    writer=Process(target=write_predictions,args=([args]))
    writer.start()
    
    #get the generator
    test_generator=get_tiledb_predict_generator(args) 
    
    #get the model
    #if calibration is to be done, get the preactivation model 
    model=get_model(args)
            
    #call the predict_on_batch_wrapper
    predict_on_batch_wrapper(args,model,test_generator)

    #drain the queue
    try:
        while not queue.empty():
            print("draining the prediction Queue")
            time.sleep(2)
    except Exception as e:
        print(e)
    
    print("joining prediction writer") 
    writer.join()
    

    
def main():
    args=predict_parse_args()
    predict(args)


if __name__=="__main__":
    main()
    
