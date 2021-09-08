from keras.utils import Sequence
import os
import signal
import psutil
import gc 
import pandas as pd
import numpy as np
import random
import math
import pysam
from ..util import *
import threading
import pickle
import pdb

def kill_child_processes(parent_pid, sig=signal.SIGTERM):
    try:
        parent = psutil.Process(parent_pid)
    except psutil.NoSuchProcess:
        return
    children = parent.children(recursive=True)
    for process in children:
        process.send_signal(sig)


def get_weights(data):
    w1=[float(data.shape[0])/sum(data.iloc[:,i]==1) for i in range(data.shape[1])]
    w0=[float(data.shape[0])/sum(data.iloc[:,i]==0) for i in range(data.shape[1])]
    return w1,w0

def open_data_file(data_path=None,tasks=None,chroms_to_use=None):
    print("running open_data_file with tasks:"+str(tasks))
    if data_path.endswith('.hdf5'):
        if tasks==None:
            data=pd.read_hdf(data_path)
        else:
            data=pd.read_hdf(data_path,columns=['CHR','START','END']+tasks)
    else:
        #treat as bed file
        if tasks==None:
            data=pd.read_csv(data_path,header=0,sep='\t')
        else:
            data=pd.read_csv(data_path,header=0,sep='\t',nrows=1)
            chrom_col=data.columns[0]
            start_col=data.columns[1]
            end_col=data.columns[2]
            data=pd.read_csv(data_path,header=0,sep='\t',usecols=[chrom_col,start_col,end_col]+tasks)
    print("loaded labels")
    print("chroms_to_use:"+str(chroms_to_use))
    print(data.head())
    try:
        data=data.set_index(['CHR','START','END'])
        print('set index to CHR, START, END')
    except:
        pass
    if chroms_to_use!=None:
        data=data[np.in1d(data.index.get_level_values(0), chroms_to_use)]
    print("filtered on chroms_to_use")
    print("data.shape:"+str(data.shape), data.columns)
    return data


#use wrappers for keras Sequence generator class to allow batch shuffling upon epoch end
class DataGenerator(Sequence):
    def __init__(self,
                 index_path,
                 input_path,
                 output_path,
                 num_inputs,
                 num_outputs,
                 ref_fasta=None,
                 batch_size=128,
                 add_revcomp=False,
                 index_tasks=None,
                 tasks=None,
                 shuffled_ref_negatives=False,
                 chroms_to_use=None,
                 get_w1_w0=False,
                 expand_dims=True,
                 upsample_thresh_list=None,
                 upsample_ratio_list=None,
                 shuffle=True,
                 return_coords=False):
        self.return_coords=return_coords
        self.expand_dims=expand_dims
        self.shuffle=shuffle
        self.batch_size=batch_size
        self.ref_fasta=ref_fasta
        self.chroms_to_use=chroms_to_use

        
        #decide if reverse complement should be used
        self.add_revcomp=add_revcomp
        if add_revcomp==True:
            self.batch_size=int(batch_size/2)

        #determine whether negative set should consist of the shuffled refs.
        # If so, split batch size in 2, as each batch will be augmented with shuffled ref negatives
        # in ratio equal to positives
        self.shuffled_ref_negatives=shuffled_ref_negatives
        if self.shuffled_ref_negatives==True:
            self.batch_size=int(self.batch_size/2)

        #get the index, input, and output files 
        self.index_tasks=index_tasks
        if tasks is None:
            tasks=[None]*num_inputs
        else:
            tasks=[i.split(',') for i in tasks] + [None]*(num_inputs-1)
        self.tasks=tasks
        print("TASKS:"+str(self.tasks))
        self.index_path=index_path
        self.input_path=input_path
        self.output_path=output_path
        self.num_inputs=num_inputs
        self.num_outputs=num_outputs

        
        self.file_to_pd=self.get_file_to_pd()        
        self.indices=self.file_to_pd[self.index_path]
        self.num_indices=self.indices.shape[0]
        print("indices:"+str(self.indices.head()))
        print("num_indices:"+str(self.num_indices))
        #handle task-specific weights -- this is a bit outdated and may be removed in the future. 
        if get_w1_w0==True:
            assert self.data is not None
            w1,w0=get_weights(self.data)
            self.w1=w1
            self.w0=w0
            
        #set variables needed for upsampling the positives
        self.upsample_thresh_list=upsample_thresh_list
        self.upsample_ratio_list=upsample_ratio_list
        #generate the upsampled threshold index subgroups
        print("creating upsampling logic for generator")
        print("self.upsample_thresh_list:"+str(self.upsample_thresh_list))
        if self.upsample_thresh_list is not None:
            self.get_upsampled_indices()
        else:
            self.indices=self.indices.index.tolist() 
            if self.shuffle == True:
                np.random.shuffle(self.indices)
        self.lock=threading.Lock() 
        print("generator initialized")
        
    def get_file_to_pd(self):
        '''
        make sure all input/output/index files are loaded only once in case there's overlap 
        generate a dictionary of file name to pandas data frame of file contents 
        '''
        file_to_df={}
        print(self.index_path)
        print(self.tasks)
        if self.tasks[0] is not None:
            file_to_df[self.index_path]=open_data_file(data_path=self.index_path,tasks=[ti[0] for ti in self.tasks if ti is not None],chroms_to_use=self.chroms_to_use)
        else:
            file_to_df[self.index_path]=open_data_file(data_path=self.index_path,tasks=self.index_tasks,chroms_to_use=self.chroms_to_use)
        print("got index_path df") 
        for i in range(self.num_inputs):
            cur_input=self.input_path[i]
            print(cur_input) 
            if cur_input=="seq":
                continue
            if cur_input in file_to_df:
                continue
            file_to_df[self.input_path[i]]=open_data_file(data_path=cur_input,tasks=self.tasks[i],chroms_to_use=self.chroms_to_use)
        print('got input')
        for i in range(self.num_outputs):
            cur_output=self.output_path[i]
            print(cur_output)
            if cur_output in file_to_df:
                print('skipped output reading')
                continue
            file_to_df[cur_output]=open_data_file(data_path=cur_output,tasks=self.tasks[i],chroms_to_use=self.chroms_to_use)
        return file_to_df
            
    def get_upsampled_indices(self):
        '''
        several levels of upsampling are handled
        self.upsample_thresh_list -- list of thresholds for upsampling the dataset 
        self.upsample_ratio_list -- fraction of batch to be generated at each threshold, length of this list = len(self.upsample_thresh_list)-1
        i.e. upsample_thresh_list=[0,0.01,1], upsample_ratio_list=[0.7,0.2] means 70% of samples should be in the range [0,0.1), 20% of samples should be 
        in the range [0.1,1), and remaining 10% of samples are in the range [1,inf)
        '''

        self.upsampled_coord_indices = {}
        self.upsampled_numerical_indices = {}
        self.batch_sizes = []

        print("upsample thresh list:"+str(self.upsample_thresh_list))
        if self.upsample_ratio_list is None:
            #all thresholds represented equally in the batch 
            self.upsample_ratio = 1 / (len(self.upsample_thresh_list) - 1)
            self.upsample_ratio_list = [self.upsample_ratio for i in range(len(self.upsample_thresh_list) - 1)]
            print("upsample ratio list: " + str(self.upsample_ratio_list))
            
        for ind in range(len(self.upsample_thresh_list)-1):
            lower_thresh_bound=self.upsample_thresh_list[ind]
            upper_thresh_bound=self.upsample_thresh_list[ind+1]
            
            #get the sub-batch that contains values within this value threshold
            sub_batch_size=int(self.batch_size*self.upsample_ratio_list[ind])
            self.batch_sizes.append(sub_batch_size)

            #get the coordinates where all values fall in the range [lower_thresh_bound, upper_thresh_bound)
            sub_batch_coords=self.indices.loc[(self.indices>=lower_thresh_bound).any(axis=1) & (self.indices < upper_thresh_bound).all(axis=1)].index
            len_sub_batch_coords=len(sub_batch_coords)
            self.upsampled_coord_indices[ind]=sub_batch_coords
            self.upsampled_numerical_indices[ind] = np.arange(len_sub_batch_coords)
            
            #number of times the current sub-set of upsampled indices should be wrapped to get to num_indices values 
            num_wraps=math.ceil(self.num_indices/len_sub_batch_coords)
            self.upsampled_numerical_indices[ind] = np.tile(self.upsampled_numerical_indices[ind], num_wraps)[0:self.num_indices]
            #shuffle the sub-set of indices, if specified 
            if self.shuffle == True:
                np.random.shuffle(self.upsampled_numerical_indices[ind])
                
        #handle the final index (i.e. unspecified upper bound)
        ind=len(self.upsample_thresh_list)-1
        lower_thresh_bound=self.upsample_thresh_list[ind]
        sub_batch_size=int(self.batch_size-sum(self.batch_sizes))
        self.batch_sizes.append(sub_batch_size)
        sub_batch_coords=self.indices.loc[(self.indices>=lower_thresh_bound).any(axis=1)].index
        len_sub_batch_coords=len(sub_batch_coords)
        self.upsampled_coord_indices[ind]=sub_batch_coords
        self.upsampled_numerical_indices[ind] = np.arange(len_sub_batch_coords)        
        #number of times the current sub-set of upsampled indices should be wrapped to get to num_indices values 
        num_wraps=math.ceil(self.num_indices/len_sub_batch_coords)
        self.upsampled_numerical_indices[ind] = np.tile(self.upsampled_numerical_indices[ind], num_wraps)[0:self.num_indices]
        #shuffle the sub-set of indices, if specified 
        if self.shuffle == True:
            np.random.shuffle(self.upsampled_numerical_indices[ind])
        del self.indices 
        return
            
        
    def __len__(self):
        return math.ceil(self.num_indices/self.batch_size)

    def get_coords(self,idx):
        if self.upsample_thresh_list is not None:
            all_bed_entries=[]
            for ind,val in enumerate(self.batch_sizes):
                batch_indices = self.upsampled_numerical_indices[ind][idx*val:(idx+1)*val]
                bed_entries = self.upsampled_coord_indices[ind][batch_indices].tolist()
                all_bed_entries+=bed_entries
        else:
            all_bed_entries=self.indices[idx*self.batch_size:(idx+1)*self.batch_size]
        return all_bed_entries
    
    def get_seq(self,coords):
        seqs=[self.ref.fetch(i[0],i[1],i[2]) for i in coords]
        return seqs
        
    def get_pd_vals(self,coords,io_index):
        try:
            return self.file_to_pd[io_index].loc[coords].values
        except:
            raise Exception("could not fetch coords:"+str(coords))
        
    def transform_seq(self,seqs):
        if self.add_revcomp==True:
            #add in the reverse-complemented sequences for training.
            seqs_rc=[revcomp(s) for s in seqs]
            seqs=seqs+seqs_rc
        if self.shuffled_ref_negatives is True:
            #generate the corresponding negative set by dinucleotide-shuffling the sequences
            seqs_shuffled=[dinuc_shuffle(s) for s in seqs]
            seqs=seqs+seqs_shuffled
        return seqs
    
    
    def transform_vals(self,vals):
        if self.add_revcomp==True:
            vals=np.concatenate((vals,vals),axis=0)
        if self.shuffled_ref_negatives is True: 
            val_shape=vals.shape
            vals=np.concatenate((vals,np.zeros(vals_shape)))
        return vals
        
    def __getitem__(self,idx):
        
        try:
            gc.unfreeze()
            self.lock.acquire() 
            #print("STARTING")
            self.ref=pysam.FastaFile(self.ref_fasta)
            #get the coordinates for the current batch
            coords=self.get_coords(idx)
            #print("GOT COORDS") 
            #get the inputs
            X=[]
            for cur_input_index in range(self.num_inputs):
                cur_input=self.input_path[cur_input_index]
                if cur_input=="seq":
                    cur_x=one_hot_encode(self.transform_seq(self.get_seq(coords)))
                    if self.expand_dims==True:
                        cur_x=np.expand_dims(cur_x,axis=1)
                else:
                    #extract values from pandas df 
                    cur_x=self.transform_vals(self.get_pd_vals(coords,cur_input))
                X.append(cur_x)
            #get the outputs
            y=[]
            for cur_output_index in range(self.num_outputs):
                cur_output=self.output_path[cur_output_index] 
                if cur_output=="seq":
                    cur_y=one_hot_encode(self.transform_seq(self.get_seq(coords)))
                    if self.expand_dims==True:
                        cur_y=np.expand_dims(cur_y,axis=1)
                else:
                    cur_y=self.transform_vals(self.get_pd_vals(coords,cur_output))
                y.append(cur_y)
            self.lock.release() 
            #return the batch as an X,y tuple
            #print("SUCCESS") 
            if self.return_coords is False: 
                return (X,y)
            else:
                return (X,y,coords)
        except Exception as e:
            print(str(e)+" from id:"+str(idx))
            kill_child_processes(os.getpid())
            raise 




    def on_epoch_end(self):
        #if upsampling is being used, shuffle the positive and negative indices
        if self.shuffle==True:
            if self.upsample_thresh_list is not None:
                for ind,val in enumerate(self.batch_sizes):
                    np.random.shuffle(self.upsampled_numerical_indices[ind])
            else:
                np.random.shuffle(self.indices)

