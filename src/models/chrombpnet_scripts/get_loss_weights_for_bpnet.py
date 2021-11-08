import tiledb
from kerasAC.tiledb_config import * 
from multiprocessing import Pool
import argparse
import os
import signal
import psutil
import numpy as np
import matplotlib.pyplot as plt
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


def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--tdb_array")
    parser.add_argument("--chroms",nargs="+")
    parser.add_argument("--ambig_attribute",default=None)
    parser.add_argument("--upsample_attribute",default=None)
    parser.add_argument("--label_attribute")
    parser.add_argument("--num_threads",type=int,default=1)
    parser.add_argument("--task",default=None)
    parser.add_argument("--task_index",default=None)
    parser.add_argument("--upsample_thresh",type=float,default=None)
    parser.add_argument("--flank",type=int,default=500)
    parser.add_argument("--sample_size",type=int,default=None)
    parser.add_argument("--outf",type=str)
    return parser.parse_args()
    

def get_chrom_index_ranges(array,chroms_to_use):
    '''
    find tdb indices corresponding to the used chromosomes 
    '''
    num_chroms=array.meta['num_chroms']
    chrom_indices=[]
    for i in range(num_chroms):
        chrom_name=array.meta['chrom_'+str(i)]
        if chrom_name in chroms_to_use:
            start_index=array.meta['offset_'+str(i)]
            end_index=start_index+array.meta['size_'+str(i)]
            chrom_indices.append((start_index,end_index))
    return chrom_indices 


def get_task_index(array,task):
    '''
    get tdb indices of user-specified tasks 
    '''
    num_tasks=array.meta['num_tasks']
    for i in range(num_tasks):
        cur_task=array.meta['task_'+str(i)]
        if cur_task==task:
            return i

def get_region_counts(inputs):
    start_index=inputs[0][0]
    end_index=inputs[0][1]
    ambig_attribute=inputs[3]
    label_attribute=inputs[2]
    task_index=inputs[4]
    upsample_thresh=inputs[5]
    upsample_attribute=inputs[6]
    flank=inputs[7]
    sample_size=inputs[8]
    counts=[]
    print("ambig_attribute:"+str(ambig_attribute))
    print("label_attribute:"+str(label_attribute))
    print("upsample_attribute:"+str(upsample_attribute))
    with tiledb.open(inputs[1],mode='r',ctx=tiledb.Ctx(get_default_config())) as array:
        if ambig_attribute is not None:
            print("starting query for ambig vals:"+str(start_index)+":"+str(end_index))
            ambig_vals=array.query(attrs=[ambig_attribute])[start_index:end_index-1,task_index][ambig_attribute]
        if upsample_attribute is not None:
            print("starting query for upsample vals:"+str(start_index)+":"+str(end_index))
            upsample_vals=array.query(attrs=[upsample_attribute])[start_index:end_index-1,task_index][upsample_attribute]
        print("starting query for label vals:"+str(start_index)+":"+str(end_index))
        label_vals=array.query(attrs=[label_attribute])[start_index:end_index-1,task_index][label_attribute]
        print("completed queries")
        if (ambig_attribute is not None) and (upsample_attribute is not None) :
            indices_for_training=np.where(np.logical_and(ambig_vals == 0, upsample_vals >= upsample_thresh))[0]
        elif (upsample_attribute is not None): 
            indices_for_training=np.where(upsample_vals >= upsample_thresh)[0]
        elif (ambig_attribute is not None):
            non_ambig_indices=np.where(ambig_vals ==0)[0]
            np.random.seed(1234)
            indices_for_training=np.random.choice(non_ambig_indices,sample_size)
        else:
            np.random.seed(1234)
            indices_for_training=np.random.choice(np.arange(label_vals.shape[0]),sample_size)
    print("got indices for region")
    for index in indices_for_training:
        if index%1000==0:
            print(str(index))
        try:
            counts.append(np.sum(label_vals[index-flank:index+flank]))
        except:
            #ran off array edge 
            continue
    return counts

def get_counts_loss_weight(tdb_path,chroms,ambig_attribute,label_attribute,upsample_thresh,upsample_attribute,flank,threads=1,task=None,task_index=None,sample_size=None):
    array=tiledb.open(tdb_path,mode='r')
    print("opened array:"+str(tdb_path) + " for reading")
    if task is not None:
        task_index=get_task_index(array,task)
    tdb_indices=get_chrom_index_ranges(array,chroms)
    pool_inputs=[]
    for entry in tdb_indices:
        pool_inputs.append([entry,tdb_path,label_attribute,ambig_attribute,task_index,upsample_thresh,upsample_attribute,flank,sample_size])
    print("got tdb indices and pool inputs")
    pool=Pool(processes=threads,initializer=init_worker)
    counts=None
    try:
        for region_counts in pool.map(get_region_counts,pool_inputs):
            if counts is None:
                counts=region_counts
            else:
                counts=np.concatenate((counts,region_counts))
        pool.close()
        pool.join()
        #summarize counts
        median_counts=np.median(counts)
        scaled_counts=median_counts/10
        #counts_threshold=np.log(np.amin(counts)/2+0.001)
        counts_threshold=np.log(np.quantile(counts, 0.1)/2+0.001)
        return scaled_counts, counts_threshold
    except KeyboardInterrupt:
        kill_child_processes(os.getpid())
        pool.terminate()
        raise
    except Exception as e:
        print(e)
        kill_child_processes(os.getpid())
        raise 

def main():
    args=parse_args()
    counts_loss_weight, counts_threshold=get_counts_loss_weight(tdb_path=args.tdb_array,
                                              chroms=args.chroms,
                                              ambig_attribute=args.ambig_attribute,
                                              label_attribute=args.label_attribute,
                                              task=args.task,
                                              task_index=args.task_index,
                                              threads=args.num_threads,
                                              upsample_thresh=args.upsample_thresh,
                                              upsample_attribute=args.upsample_attribute,
                                              flank=args.flank,
                                              sample_size=args.sample_size)
    #print("counts_loss_weight:"+str(counts_loss_weight))
    f = open(args.outf, "w")
    f.write(str(counts_loss_weight))
    f.write("\n")
    f.write(str(counts_threshold))
    f.close()

if __name__=="__main__":
    main()
    
