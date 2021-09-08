from tensorflow.keras.utils import Sequence
import itertools
import os
import signal
import psutil
import pandas as pd
import numpy as np
from scipy.special import logit 
import random
import threading 
from random import shuffle 
import math
from math import ceil, floor
import pysam
from ..util import *
from ..tiledb_config import * 
import tiledb
import pdb
from ..s3_sync import * 
from collections import OrderedDict
import gc
import pdb         
    
import numpy as np
import tensorflow as tf
import random as rn
import os
from scipy.special import softmax


os.environ['PYTHONHASHSEED'] = '0'
np.random.seed(0)
rn.seed(0)
tf.random.set_seed(0)

def get_upsampled_indices_chrom(inputs):
    region_start=inputs[0]
    region_end=inputs[1]
    tdb_array_name=inputs[2]
    tdb_ambig_attribute=inputs[3]
    tdb_partition_attribute_for_upsample=inputs[4]
    dataset_indices=inputs[5]
    tdb_partition_thresh_for_upsample=inputs[6]
    print("starting getting indices to upsample in range:"+str(region_start)+"-"+str(region_end))
    with tiledb.open(tdb_array_name,'r',ctx=tiledb.Ctx(get_default_config())) as tdb_array:
        if tdb_ambig_attribute is not None:
            attr_vals=tdb_array.query(attrs=[tdb_ambig_attribute,tdb_partition_attribute_for_upsample]).multi_index[region_start:region_end-1,dataset_indices]
            ambig_attr_vals=np.sum(attr_vals[tdb_ambig_attribute],axis=1)
        else:
            attr_vals=tdb_array.query(attrs=[tdb_partition_attribute_for_upsample]).multi_index[region_start:region_end-1,dataset_indices]        
        upsample_vals=np.sum(attr_vals[tdb_partition_attribute_for_upsample],axis=1)
    if tdb_ambig_attribute is not None:
        cur_upsampled_indices=region_start+np.argwhere((upsample_vals>=tdb_partition_thresh_for_upsample) & ( ambig_attr_vals==0))
    else: 
        cur_upsampled_indices=region_start+np.argwhere(upsample_vals>=tdb_partition_thresh_for_upsample)
    print("finished indices to upsample in range:"+str(region_start)+"-"+str(region_end))
    return cur_upsampled_indices

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



class TiledbGenerator(Sequence):
    def __init__(self,
                 ref_fasta,
                 batch_size,
                 tdb_array,
                 tdb_partition_attribute_for_upsample,
                 tdb_partition_thresh_for_upsample,
                 tdb_partition_datasets_for_upsample,
                 tdb_input_source_attribute,
                 tdb_input_flank,
                 upsample_ratio,
                 tdb_output_source_attribute,
                 tdb_output_flank,
                 num_inputs,
                 num_outputs,
                 tdb_input_min=None,
                 tdb_input_max=None,
                 tdb_output_min=None,
                 tdb_output_max=None,
                 tdb_input_datasets=None,
                 tdb_output_datasets=None,
                 tdb_input_aggregation=None,
                 tdb_input_transformation=None,
                 tdb_output_aggregation=None,
                 tdb_output_transformation=None,
                 tdb_ambig_attribute=None,
                 chroms=None,
                 chrom_sizes=None,
                 shuffle_epoch_start=True,
                 shuffle_epoch_end=True,
                 pseudocount=0.001,
                 add_revcomp=False,
                 return_coords=False,
                 tdb_config=None,
                 tdb_ctx=None,
                 num_threads=1,
                 bed_regions=None,
                 bed_regions_center=None,
                 bed_regions_jitter=1):
        '''
        tdb_partition_attribute_for_upsample -- attribute in tiledb array used for determining which bases to upsample (usu. 'idr_peak') 
        tdb_partition_thresh_for_upsample -- threshold for determinining samples to upsample (generally 1) 
        tdb_input_aggregation/ tdb_output_aggregation -- one of 'average','max','binary_max','sum',None
        '''
        self.num_threads=num_threads
        self.shuffle_epoch_start=shuffle_epoch_start
        self.shuffle_epoch_end=shuffle_epoch_end

        #get local copy of s3 reference sequence
        if ref_fasta.startswith('s3://'):
            self.ref_fasta=download_s3_file(ref_fasta)
            fai=download_s3_file(ref_fasta+'.fai')
        else: 
            self.ref_fasta=ref_fasta

        self.batch_size=batch_size

        self.add_revcomp=add_revcomp
        if self.add_revcomp==True:
            self.batch_size=int(math.floor(self.batch_size/2))

        #create tiledb configuration parameters (these have been found optimal for most use cases, but should set in a separate config file in the future)
        if tdb_config is not None:
            self.config=tdb_config
        else:
            self.config=get_default_config()
        if tdb_ctx is not None:
            self.ctx=tdb_ctx
        else:
            self.ctx=tiledb.Ctx(self.config)
            
        print("opening:"+tdb_array+" for reading...")
        self.tdb_array_name=tdb_array
        self.tdb_array=tiledb.open(tdb_array,mode='r',ctx=self.ctx)
        
        print("success!")

        #identify chromosome information
        if chroms is not None:
            self.chroms_to_use=chroms
        else:
            if chrom_sizes.startswith("s3://"):
                self.chroms_to_use=[i.split()[0] for i in read_s3_file_contents(chrom_sizes).strip().split('\n')]
            else:
                self.chroms_to_use=[i.split()[0] for i in open(chrom_sizes,'r').read().strip().split('\n')]
        #find the tdb indices that correspond to the chroms to be used 
        self.get_chrom_index_ranges(self.chroms_to_use)
        #print("self.weighted_chrom_indices"+str(self.weighted_chrom_indices))
        print("got indices for used chroms")


        #get indices of input datasts to be used in training
        self.input_dataset_indices=self.get_dataset_indices(tdb_input_datasets)
        print("identified input dataset indices:"+str(self.input_dataset_indices))

        #get indices for outputs to be used in training
        self.output_dataset_indices=self.get_dataset_indices(tdb_output_datasets)
        print("identified output dataset indices:"+str(self.output_dataset_indices))

        self.tdb_ambig_attribute=tdb_ambig_attribute

        #store input params
        self.num_inputs=num_inputs
        self.tdb_input_source_attribute=[i.split(',') for i in tdb_input_source_attribute]
        self.tdb_input_flank=[[int(j) for j in  i.split(',')] for i in tdb_input_flank]
        self.tdb_input_aggregation=[[str(j) for j in  i.split(',')] for i in tdb_input_aggregation]
        self.tdb_input_transformation=[[str(j) for j in  i.split(',')] for i in tdb_input_transformation]

        #store output params
        self.num_outputs=num_outputs
        self.tdb_output_source_attribute=[i.split(',') for i in tdb_output_source_attribute]
        self.tdb_output_flank=[[int(j) for j in  i.split(',')] for i in tdb_output_flank]
        self.tdb_output_aggregation=[[str(j) for j in  i.split(',')] for i in tdb_output_aggregation]
        self.tdb_output_transformation=[[str(j) for j in  i.split(',')] for i in tdb_output_transformation]


        #identify min/max values
        self.tdb_input_min=transform_data_type_min(tdb_input_min,self.num_inputs)
        self.tdb_input_max=transform_data_type_max(tdb_input_max,self.num_inputs)
        self.tdb_output_min=transform_data_type_min(tdb_output_min,self.num_outputs)
        self.tdb_output_max=transform_data_type_max(tdb_output_max,self.num_outputs)
                
        #identify upsampled genome indices for model training
        self.tdb_partition_attribute_for_upsample=tdb_partition_attribute_for_upsample
        self.tdb_partition_thresh_for_upsample=tdb_partition_thresh_for_upsample
        self.tdb_partition_datasets_for_upsample=tdb_partition_datasets_for_upsample,

    
        #handle the option of training/predicting on pre-specified bed regions
        if bed_regions is not None:
            if type(bed_regions)==str: 
                self.bed_regions=pd.read_csv(bed_regions,header=None,sep='\t')
            else:
                self.bed_regions=bed_regions
            self.bed_regions_center=bed_regions_center
            if self.bed_regions_center=="random":
                self.bed_regions_jitter=bed_regions_jitter
            print("loaded bed regions")
            #get mapping of bed region to tdb index
            self.map_regions_to_tdb_index() 
        else:
            self.bed_regions=None 
            self.coord=None
            if (upsample_ratio is not None) and (upsample_ratio > 0):
                assert type(upsample_ratio)==float
            self.upsample_ratio=upsample_ratio
            if (self.upsample_ratio is not None) and (upsample_ratio >0):
                #get indices for dataset used for upsamples
                self.partition_thresh_dataset_indices=list(set(itertools.chain.from_iterable(self.get_dataset_indices(tdb_partition_datasets_for_upsample))))
                print("identified upsampling dataset indices:"+str(self.partition_thresh_dataset_indices))
                self.get_upsampled_indices()
                self.upsampled_batch_size=math.ceil(self.upsample_ratio*self.batch_size)
            else:
                self.upsampled_batch_size=0
                self.upsampled_indices_len=0
                self.upsampled_indices=[]
            
            self.non_upsampled_batch_size=self.batch_size-self.upsampled_batch_size        

        self.pseudocount=pseudocount
        self.return_coords=return_coords
            
        print('created generator')
        
    def map_regions_to_tdb_index(self):
        self.coord=[]
        self.tdb_indices=[]
        for index,row in self.bed_regions.iterrows():
            chrom=row[0]
            start=row[1]
            end=row[2]
            if self.bed_regions_center =="summit":
                summit=row[9]
                pos=start+summit
                tdb_index=self.chrom_to_indices[chrom][0]+pos
                self.coord.append([chrom,pos])
                self.tdb_indices.append(tdb_index) 
            elif self.bed_regions_center == "center":
                pos=int(round((start+end)/2))
                tdb_index=self.chrom_to_indices[chrom][0]+pos
                self.coord.append([chrom,pos])
                self.tdb_indices.append(tdb_index)                 
            elif self.bed_regions_center == "edges":
                ## this mode is introuduced to test at edges on CHIP-seq
                ## assuming the same outputlength in  multi task output
                pos=start+row[9]-self.tdb_output_flank[0]
                tdb_index=self.chrom_to_indices[chrom][0]+pos
                self.coord.append([chrom,pos])
                self.tdb_indices.append(tdb_index)                 
                pos=start+row[9]+self.tdb_output_flank[0]
                tdb_index=self.chrom_to_indices[chrom][0]+pos
                self.coord.append([chrom,pos])
                self.tdb_indices.append(tdb_index)                 
            else:
                assert self.bed_regions_center =="random"
                #select n=bed_regions_jitter bases from each peak to center the training/validation interval
                for jitter_index in range(self.bed_regions_jitter): 
                    pos=random.randint(start,end) 
                    tdb_index=self.chrom_to_indices[chrom][0]+pos
                    self.coord.append([chrom,pos])
                    self.tdb_indices.append(tdb_index)
        #shuffle the jittered bed regions, preserving correspondence of self.tdb_indices & self.coord
        temp = list(zip(self.coord, self.tdb_indices)) 
        if self.shuffle_epoch_start==True:
            random.shuffle(temp) 
        self.coord, self.tdb_indices = zip(*temp)
        
    def get_chrom_index_ranges(self,chroms_to_use):
        '''
        find tdb indices corresponding to the used chromosomes 
        '''
        num_chroms=self.tdb_array.meta['num_chroms']
        self.dataset_indices=[i for i in range(self.tdb_array.meta['num_tasks'])]
        chrom_indices=[]
        chrom_sizes=[]
        chroms=[]
        num_indices=0
        for i in range(num_chroms):
            chrom_name=self.tdb_array.meta['chrom_'+str(i)]
            if chrom_name in chroms_to_use:
                chroms.append(chrom_name)
                start_index=self.tdb_array.meta['offset_'+str(i)]
                end_index=start_index+self.tdb_array.meta['size_'+str(i)]
                num_indices+=(end_index-start_index)
                chrom_indices.append((start_index,end_index))
                chrom_sizes.append(self.tdb_array.meta['size_'+str(i)])
        min_chrom_size=min(chrom_sizes)
        scaled_chrom_sizes=[round(i/min_chrom_size) for i in chrom_sizes]
        weighted_chrom_sizes=[]
        for i in range(len(chrom_sizes)):
            cur_weight=scaled_chrom_sizes[i]
            cur_range=[chrom_indices[i]]
            weighted_chrom_sizes=weighted_chrom_sizes+cur_weight*cur_range
        self.chrom_indices=chrom_indices
        self.weighted_chrom_indices=weighted_chrom_sizes
        self.num_indices=num_indices
        self.chroms_to_use=chroms
        self.chrom_to_indices={}
        for i in range(len(self.chroms_to_use)):
            cur_chrom=self.chroms_to_use[i]
            cur_indices=self.chrom_indices[i]
            self.chrom_to_indices[cur_chrom]=cur_indices 
        return
    
    def get_dataset_indices(self,dataset_names):
        '''
        get tdb indices of user-specified tasks. 
        returns a list of lists -- inner list refers to tasks, outer list refers to either inputs or outputs. 
        '''
        num_datasets=self.tdb_array.meta['num_tasks'] #tdb array still uses num_tasks in metadata, this should eventually become num_datasets
        dataset_indices=[]
        for io_index in range(len(dataset_names)):
            dataset_indices.append([])
            datasets_by_task=dataset_names[io_index].split(',')
            for task_index in range(len(datasets_by_task)):
                cur_io_cur_task_dataset=datasets_by_task[task_index]
                for i in range(num_datasets):
                    tdb_dataset=self.tdb_array.meta['task_'+str(i)]
                    if tdb_dataset == cur_io_cur_task_dataset:
                        dataset_indices[io_index].append(i)
        assert(len(dataset_indices)>0)
        return dataset_indices
    
    def get_nonupsample_batch_indices(self):
        '''
        randomly select n positions from the genome 
        '''
        #get current chromosome
        cur_interval=random.sample(self.weighted_chrom_indices,1)[0]
        #sample random indices from the current chromosome 
        cur_batch=random.sample(range(cur_interval[0],cur_interval[1]),self.non_upsampled_batch_size)
        return cur_batch
                                                                                                    
    

    def get_upsampled_indices(self):
        from multiprocessing import Pool
        print("num_threads:"+str(self.num_threads))
        pool=Pool(processes=self.num_threads,initializer=init_worker)
        pool_inputs=[] 
        for region in self.chrom_indices:
            region_start=region[0]
            region_end=region[1]
            pool_inputs.append((region_start,region_end,self.tdb_array_name,self.tdb_ambig_attribute,self.tdb_partition_attribute_for_upsample,self.partition_thresh_dataset_indices,self.tdb_partition_thresh_for_upsample))
        upsampled_indices=None
        try:
            for region_upsampled_indices in pool.map(get_upsampled_indices_chrom,pool_inputs):
                if upsampled_indices is None:
                    upsampled_indices=np.squeeze(region_upsampled_indices)
                else:
                    upsampled_indices=np.concatenate((upsampled_indices,np.squeeze(region_upsampled_indices)))
        except KeyboardInterrupt:
            kill_child_processes(os.getpid())
            pool.terminate()
            raise
        except Exception as e:
            print(e)
            kill_child_processes(os.getpid())
            raise 
        pool.close()
        pool.join()
        print('closed upsampling pool') 
        print("made upsampled index data frame")
        self.upsampled_indices=upsampled_indices
        if self.shuffle_epoch_start==True:
            #shuffle rows & reset index
            print("shuffling upsampled dataframes prior to start of training")
            np.random.shuffle(self.upsampled_indices)
        self.upsampled_indices_len=len(self.upsampled_indices)
        print("finished upsampling")
        return

    
        
    def __len__(self):
        #we have an explict set of regions
        if self.bed_regions is not None:
            return int(ceil(len(self.tdb_indices)/self.batch_size))
        #we are only training on peak regions
        elif (self.upsample_ratio is not None) and (self.upsample_ratio==1):
            return int(ceil(self.upsampled_indices_len/self.upsampled_batch_size))
        else:
        #training on peak and non-peak regions 
            return int(ceil(self.num_indices/self.batch_size))
    

    def __getitem__(self,idx):
        gc.unfreeze()
        self.ref=pysam.FastaFile(self.ref_fasta)
        
        #get the coordinates for the current batch
        tdb_batch_indices=self.get_tdb_indices_for_batch(idx) #coords is a df with 'chrom' and 'pos' columns.
        coords=None
        if self.return_coords is True:
            #get the chromosome coordinates that correspond to indices
            coords=self.get_coords(tdb_batch_indices,idx)
        #get the inputs 
        X=[]
        #iterate through the list of model inputs 
        for cur_input_index in range(self.num_inputs):
            cur_input=self.tdb_input_source_attribute[cur_input_index]
            cur_x=None
            #iterate through the stacked channels of the current input
            for cur_input_channel_index in range(len(cur_input)):
                cur_input_channel=cur_input[cur_input_channel_index]
                if cur_input_channel=="seq":                
                    #get the one-hot encoded sequence
                    if coords is None:
                        coords=self.get_coords(tdb_batch_indices,idx)
                    cur_seq=one_hot_encode(self.get_seq(coords,self.tdb_input_flank[cur_input_index][cur_input_channel_index]))
                    if cur_x is None: 
                        cur_x=cur_seq
                    else:
                        cur_x=np.concatenate((cur_x,cur_seq),axis=-1)
                else:
                    #extract values from tdb
                    #print(tdb_batch_indices)
                    #print(cur_input_channel,self.tdb_input_flank,cur_input_index,cur_input_channel_index,self.input_dataset_indices)
                    cur_vals=self.get_tdb_vals(tdb_batch_indices,cur_input_channel,self.tdb_input_flank[cur_input_index][cur_input_channel_index],self.input_dataset_indices[cur_input_index][cur_input_channel_index])
                    aggregate_vals=self.aggregate_vals(cur_vals,self.tdb_input_aggregation[cur_input_index][cur_input_channel_index])
                    transformed_vals=self.transform_vals(aggregate_vals,self.tdb_input_transformation[cur_input_index][cur_input_channel_index])
                    if cur_x is None: 
                        cur_x=transformed_vals
                    else:
                        cur_x=np.concatenate((cur_x,transformed_vals),axis=-1)
            
            #perform reverse complementation, if specified
            if self.add_revcomp is True:
                cur_x=np.concatenate((cur_x,np.flip(cur_x)),axis=0)
            X.append(cur_x)
                 
        #get the outputs 
        y=[]
        for cur_output_index in range(self.num_outputs):
            cur_y=None
            cur_output=self.tdb_output_source_attribute[cur_output_index]
            for cur_output_channel_index in range(len(cur_output)):
                cur_output_channel=cur_output[cur_output_channel_index]
                if cur_output_channel=="seq":
                    #get the one-hot encoded sequence
                    if coords is None:
                        coords=get_coords(tdb_batch_indices,idx)
                    cur_seq=one_hot_encode(self.get_seq(coords,self.tdb_output_flank[cur_output_index][cur_output_channel_index]))
                    if cur_y is None:
                        cur_y=cur_seq
                    else:
                        cur_y=np.concatenate((cur_y,cur_seq),axis=-1)
                else:
                    #extract values from tdb
                    cur_vals=self.get_tdb_vals(tdb_batch_indices,cur_output_channel,self.tdb_output_flank[cur_output_index][cur_output_channel_index],self.output_dataset_indices[cur_output_index][cur_output_channel_index])
                    aggregate_vals=self.aggregate_vals(cur_vals,self.tdb_output_aggregation[cur_output_index][cur_output_channel_index])
                    transformed_vals=self.transform_vals(aggregate_vals,self.tdb_output_transformation[cur_output_index][cur_output_channel_index])
                    if cur_y is None:
                        cur_y=transformed_vals
                    else:
                        cur_y=np.concatenate((cur_y,transformed_vals),axis=-1)
            if self.add_revcomp is True:
                cur_y=np.concatenate((cur_y,np.flip(cur_y)),axis=0)
            y.append(cur_y)
            
        if self.return_coords is True:
            coords_updated=[]
            if self.add_revcomp==True:
                for i in coords:
                    coords_updated.append(i+['p'])
                for i in coords:
                    coords_updated.append(i+['r'])
            else:
                for i in coords:
                    coords_updated.append(i+['.'])
            coords= np.string_(coords_updated)
            coords = coords.astype('S256')            
            filtered_X,filtered_y,filtered_coords=self.remove_data_out_of_range(X,y,coords)
        else:
            filtered_X,filtered_y,filtered_coords=self.remove_data_out_of_range(X,y)
        if filtered_X[0].size==0:
            #empty!
            try:
                return self.__getitem__(idx+1)
            except:
                #we are at the last index, wrap around 
                return self.__getitem__(0)
        if self.return_coords is True:
            #print(str(filtered_coords))
            return (filtered_X,filtered_y,filtered_coords)
        else:
            return (filtered_X,filtered_y)
        
    def remove_data_out_of_range(self,X,y,coords=None):
        bad_indices=[]
        for i in range(len(self.tdb_input_min)):
            out_of_range=[z[0] for z in np.argwhere(X[i]<self.tdb_input_min[i]).tolist() if len(z)>0 ]
            bad_indices+=out_of_range
            out_of_range=[z[0] for z in np.argwhere(X[i]>self.tdb_input_max[i]).tolist() if len(z)>0 ]
            bad_indices+=out_of_range
        for i in range(len(self.tdb_output_min)):
            out_of_range=[z[0] for z in np.argwhere(y[i]<self.tdb_output_min[i]).tolist() if len(z)>0]
            bad_indices+=out_of_range
            out_of_range=[z[0] for z in np.argwhere(y[i]>self.tdb_output_max[i]).tolist() if len(z)>0]
            bad_indices+=out_of_range
        bad_indices=list(set(bad_indices))
        X=[np.delete(i,bad_indices,0) for i in X]
        y=[np.delete(i,bad_indices,0) for i in y]
        if coords is not None:
            coords=np.delete(coords,bad_indices,0)
        return X,y,coords
        
    def get_coords(self,tdb_batch_indices,idx):
        #return list of (chrom,pos) for each index in batch
        
        #if we are using bed regions supplied by user, the coords have already been pre-computed 
        if self.coord is not None:
            return self.coord[idx*self.batch_size:(idx+1)*self.batch_size]
        coords=[]
        for cur_batch_index in tdb_batch_indices:
            for chrom_index in range(len(self.chrom_indices)):
                cur_chrom_start_index=self.chrom_indices[chrom_index][0]
                cur_chrom_end_index=self.chrom_indices[chrom_index][1]
                if (cur_batch_index >=cur_chrom_start_index) and (cur_batch_index<cur_chrom_end_index):
                    coords.append([self.chroms_to_use[chrom_index],cur_batch_index-cur_chrom_start_index])
                    break
        return coords
    
    def get_tdb_indices_for_batch(self,idx):
        if self.bed_regions is not None:
            return self.tdb_indices[idx*self.batch_size:(idx+1)*self.batch_size]
        
        upsampled_batch_indices=None
        non_upsampled_batch_indices=None
        if self.upsampled_batch_size > 0:
            #might need to wrap to get the upsampled index length
            upsampled_batch_start=int(idx*self.upsampled_batch_size % self.upsampled_indices_len)
            upsampled_batch_end=upsampled_batch_start+self.upsampled_batch_size
            while upsampled_batch_end > self.upsampled_indices_len:
                if upsampled_batch_indices is None:
                    upsampled_batch_indices=self.upsampled_indices[upsampled_batch_start:self.upsampled_indices_len]
                else:
                    upsampled_batch_indices=np.concatenate((upsampled_batch_indices,self.upsampled_indices[upsampled_batch_start:self.upsampled_indices_len]))
                upsampled_batch_start=0
                upsampled_batch_end=upsampled_batch_end-self.upsampled_indices_len
            if upsampled_batch_indices is None:
                upsampled_batch_indices=self.upsampled_indices[upsampled_batch_start:upsampled_batch_end]
            else: 
                upsampled_batch_indices=np.concatenate((upsampled_batch_indices,self.upsampled_indices[upsampled_batch_start:upsampled_batch_end]))
            
        if self.non_upsampled_batch_size > 0:
            #select random indices from genome
            non_upsampled_batch_indices=self.get_nonupsample_batch_indices()
        if (upsampled_batch_indices is not None) and (non_upsampled_batch_indices is not None):
            tdb_batch_indices=np.concatenate((upsampled_batch_indices,non_upsampled_batch_indices))
        elif upsampled_batch_indices is not None:
            tdb_batch_indices=upsampled_batch_indices
        elif non_upsampled_batch_indices is not None:
            tdb_batch_indices=non_upsampled_batch_indices
        else:
            raise Exception("both upsampled_batch_indices and non_upsampled_batch_indices appear to be none")
        return tdb_batch_indices
    
     
    def get_seq(self,coords,flank):
        seqs=[]
        for coord in coords:
            chrom=coord[0]
            start_pos=coord[1]-flank
            end_pos=coord[1]+flank 
            try:
                seq=self.ref.fetch(chrom,start_pos,end_pos)
                if len(seq)<2*flank:
                    delta=2*flank-len(seq)
                    seq=seq+"N"*delta
            except:
                seq="N"*2*flank
            seqs.append(seq) 
        return seqs

    '''
    def get_tdb_vals(self,tdb_batch_indices,attribute,flank,dataset_index):
        num_entries=len(tdb_batch_indices)
        pdb.set_trace()
        #prepopulate the values array with nans
        vals=np.full((num_entries,2*flank,1),np.nan)
        #iterate through entries
        for val_index in range(num_entries):
            vals[val_index,:,:]=self.tdb_array.query(attrs=[attribute]).multi_index[tdb_batch_indices[val_index]-flank:tdb_batch_indices[val_index]+flank-1,dataset_index][attribute]
        return vals
    '''

    def get_tdb_vals(self,tdb_batch_indices,attribute,flank,dataset_index):
        flattened_batch_indices=[slice(i-flank,i+flank-1) for i in tdb_batch_indices]
        vals=self.tdb_array.query(attrs=[attribute]).multi_index[flattened_batch_indices,dataset_index][attribute]
        vals=np.reshape(vals,(len(tdb_batch_indices),-1))
        vals=np.expand_dims(vals,axis=-1)
        return vals

    
    def transform_vals(self,vals,transformer):
        if transformer is None:
            return vals 
        if transformer == 'None':
            return vals
        elif transformer == 'asinh':
            return np.arcsinh(vals)
        elif transformer == 'round':
            return np.around(vals)
        elif transformer == 'softmax':
            return softmax(vals+self.pseudocount, axis=1)
        elif transformer == 'log10':
            return np.log10(vals+self.pseudocount)
        elif transformer == 'log':
            return np.log(vals+self.pseudocount)
        elif transformer == 'counts_to_logit':
            try:
                vals=vals/np.expand_dims(vals.sum(axis=1),axis=1) #transform to probability space, axis 0 = batch, axis 1 = genome pos, axis 2 = task
            except ZeroDivisionError:
                vals=vals/(np.expand_dims(vals.sum(axis=1),axis=1)+self.pseudocount) #transform to probability space, axis 0 = batch, axis 1 = genome pos, axis 2 = task
            vals=logit(vals+self.pseudocount)
            return vals 
        else:
            raise Exception("transform_vals argument must be one of None, asinh, log10, log; you provided:"+transformer) 
    
    def aggregate_vals(self,vals,aggregator):
        if aggregator is None:
            return vals
        if aggregator == 'None':
            return vals
        elif aggregator == 'average':
            return np.mean(vals,axis=1)
        elif aggregator == 'max':
            return np.max(vals,axis=1)
        elif aggregator == 'binary_max':
            #get the max in the interval, but cap it at one or 0
            raw_max=np.max(vals,axis=1)
            raw_max[raw_max>1]=1
            raw_max[raw_max<0]=0
            return raw_max
        elif aggregator == 'sum':
            return np.sum(vals,axis=1) 
        else:
            raise Exception("aggregate_vals argument must be one of None, average, max, sum; you provided:"+aggregator)

    
    def on_epoch_end(self):
        if self.shuffle_epoch_end==True:
            if self.bed_regions is not None:
                temp = list(zip(self.coord, self.tdb_indices))
                random.shuffle(temp)
                self.coord, self.tdb_indices = zip(*temp)
            else:
                #print("WARNING: SHUFFLING ON EPOCH END MAYBE SLOW:"+str(self.upsampled_indices.shape))
                #self.upsampled_indices=self.upsampled_indices.sample(frac=1)
                np.random.shuffle(self.upsampled_indices)

