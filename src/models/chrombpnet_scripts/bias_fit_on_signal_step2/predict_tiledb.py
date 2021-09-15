
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
from kerasAC.calibrate import * 
from kerasAC.generators.basic_generator import *
from kerasAC.generators.tiledb_predict_generator import *
from kerasAC.tiledb_config import *
from kerasAC.s3_sync import *
from get_model import *
from kerasAC.splits import *
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
from keras.utils import multi_gpu_model
from kerasAC.custom_losses import *
from abstention.calibration import PlattScaling, IsotonicRegression 
import random
import pdb 

def parse_args():
    parser=argparse.ArgumentParser(description='Provide model files  & a dataset, get model predictions')
    parser.add_argument('--ref_fasta')

    tiledbgroup=parser.add_argument_group("tiledb")
    tiledbgroup.add_argument("--tdb_array",help="name of tdb array to use")
    
    tiledbgroup.add_argument("--tdb_output_source_attribute",nargs="+",default=None,help="tiledb attribute for use in label generation i.e. fc_bigwig")
    tiledbgroup.add_argument("--tdb_output_datasets",nargs="+",default=None,help="dataset column from db_ingest; comma separated across channels; space separated across outputs")    
    tiledbgroup.add_argument("--tdb_output_min",nargs="*", default=None)
    tiledbgroup.add_argument("--tdb_output_max",nargs="*", default=None)        
    tiledbgroup.add_argument("--tdb_output_flank",nargs="+",help="flank around bin center to use in generating outputs")
    tiledbgroup.add_argument("--tdb_output_aggregation",nargs="+",default=None,help="method for output aggreagtion; one of None, 'avg','max'")
    tiledbgroup.add_argument("--tdb_output_transformation",nargs="+",default=None,help="method for output transformation; one of None, 'log','log10','asinh'")

    tiledbgroup.add_argument("--tdb_input_source_attribute",nargs="+",help="attribute to use for generating model input, or 'seq' for one-hot-encoded sequence")
    tiledbgroup.add_argument("--tdb_input_min",nargs="*", default=None)
    tiledbgroup.add_argument("--tdb_input_max",nargs="*", default=None)    
    tiledbgroup.add_argument("--tdb_input_datasets",nargs="+",default=None,help="dataset column from db_ingest; comma separated across channels; space separated across inputs")
    tiledbgroup.add_argument("--tdb_input_flank",nargs="+",help="length of sequence around bin center to use for input")
    tiledbgroup.add_argument("--tdb_input_aggregation",nargs="+",default=None,help="method for input aggregation; one of 'None','avg','max'")
    tiledbgroup.add_argument("--tdb_input_transformation",nargs="+",default=None,help="method for input transformation; one of None, 'log','log10','asinh'")
    tiledbgroup.add_argument("--tdb_transformation_pseudocount",type=float,default=1)
    
    tiledbgroup.add_argument("--tdb_partition_attribute_for_upsample",default="idr_peak",help="tiledb attribute to use for upsampling, i.e. idr_peak")
    tiledbgroup.add_argument("--tdb_partition_thresh_for_upsample",type=float,default=1,help="values >= partition_thresh_for_upsample within the partition_attribute_for_upsample will be upsampled during training")
    tiledbgroup.add_argument("--tdb_partition_datasets_for_upsample",nargs="+")
    
    tiledbgroup.add_argument("--upsample_ratio_list_predict",type=float,nargs="*")
    tiledbgroup.add_argument("--tdb_ambig_attribute",default=None,help="attribute indicating ambiguous regions to not train on")
    tiledbgroup.add_argument("--chrom_sizes",default=None,help="chromsizes file for use with tiledb generator")
    tiledbgroup.add_argument("--tiledb_stride",type=int,default=1)
    tiledbgroup.add_argument("--upsample_threads",type=int,default=1)
    tiledbgroup.add_argument("--revcomp",action="store_true")
    input_filtering_params=parser.add_argument_group("input_filtering_params")    
    input_filtering_params.add_argument('--predict_chroms',nargs="*",default=None)
    input_filtering_params.add_argument("--genome",default=None)
    input_filtering_params.add_argument("--fold",type=int,default=None)
    input_filtering_params.add_argument("--bed_regions",default=None) 
    input_filtering_params.add_argument('--bed_regions_center', choices=['random','center','summit'],default=None)
    input_filtering_params.add_argument("--datasets",nargs="*",default=None)
    input_filtering_params.add_argument("--dataset_indices",nargs="*",default=None)
    
    output_params=parser.add_argument_group("output_params")
    output_params.add_argument('--predictions_and_labels_hdf5',help='name of hdf5 to save predictions',default=None)
    calibration_params=parser.add_argument_group("calibration_params")
    calibration_params.add_argument("--calibrate_classification",action="store_true",default=False)
    calibration_params.add_argument("--calibrate_regression",action="store_true",default=False)        
    
    weight_params=parser.add_argument_group("weight_params")
    weight_params.add_argument('--w1',nargs="*",type=float)
    weight_params.add_argument('--w0',nargs="*",type=float)
    weight_params.add_argument("--w1_w0_file",default=None)


    model_params=parser.add_argument_group("model_params")
    model_params.add_argument('--load_model_hdf5',help='hdf5 file that stores the model')
    model_params.add_argument('--weights',help='weights file for the model')
    model_params.add_argument('--yaml',help='yaml file for the model')
    model_params.add_argument('--json',help='json file for the model')
    model_params.add_argument('--functional',default=False,help='use this flag if your model is a functional model',action="store_true")
    model_params.add_argument('--squeeze_input_for_gru',action='store_true')
    model_params.add_argument("--num_inputs",type=int)
    model_params.add_argument("--num_outputs",type=int)
    model_params.add_argument("--num_gpus",type=int,default=1)
    
    snp_params=parser.add_argument_group("snp_params")
    snp_params.add_argument('--background_freqs',default=None)
    snp_params.add_argument('--flank',default=500,type=int)
    snp_params.add_argument('--mask',default=10,type=int)
    snp_params.add_argument('--ref_col',type=int,default=None)
    snp_params.add_argument('--alt_col',type=int,default=None)

    parser.add_argument('--batch_size',type=int,help='batch size to use to make model predictions',default=50)
    return parser.parse_args()

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
    perform_calibration=args.calibrate_classification or args.calibrate_regression
    if perform_calibration==True:
        if args.calibrate_classification==True:
            print("getting logits")
            model=Model(inputs=model.input,
                               outputs=model.layers[-2].output)
        elif args.calibrate_regression==True:
            print("getting pre-relu outputs (preacts)")
            model=Model(inputs=model.input,
                        outputs=model.layers[-1].output)
            
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
    
    #sync files to s3 if needed
    if args.predictions_and_labels_hdf5.startswith('s3://'):
        #use a local version of the file and upload to s3 when finished
        out_predictions_prefix=args.predictions_and_labels_hdf5+".predictions"
        #upload outputs for all tasks
        import glob
        to_upload=[]
        for f in glob.glob(out_predictions_prefix+"*"):
            s3_path='/'.join(out_predictions_prefix.split('/')[0:-1]+[f.split('/')[-1]])
            print(f+'-->'+s3_path)
            upload_s3_file(s3_path,f)
        
    #perform calibration, if specified
    if perform_calibration is True:
        print("calibrating")
        calibrate(args)

    #clean up any s3 artifacts:
    run_cleanup()
    
def main():
    args=parse_args()
    predict(args)


if __name__=="__main__":
    main()
    
