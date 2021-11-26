from tensorflow.keras.utils import Sequence
from utils.one_hot import dna_to_one_hot as one_hot_encode
from scipy.special import softmax
from scipy.special import logit 
from utils.tiledb_generator_utils import transform_data_type_min, transform_data_type_max, get_default_config, tdb_config_params
from math import ceil, floor
import numpy as np
import tensorflow as tf
import pandas as pd
import numpy as np
import tiledb
import pdb
import gc
import pdb
import pysam
import math
import itertools
import os
import signal
import psutil
import random

# add remove outliers and add remove blacklist regions

os.environ['PYTHONHASHSEED'] = '0'

class TiledbGenerator(Sequence):
    def __init__(self,
                 peak_regions,
                 nonpeak_regions,
                 genome_fasta,
                 batch_size,
                 input_len,
                 output_len,
                 max_jitter,
                 negative_sampling_ratio,
                 cts_sum_min_thresh,
                 cts_sum_max_thresh,
                 tdb_array,
                 tdb_dataset,
                 seed,
                 add_revcomp,
                 return_coords
                 ):

        np.random.seed(seed)
        random.seed(seed)
        tf.random.set_seed(seed)  

        self.peak_regions=peak_regions
        self.nonpeak_regions=nonpeak_regions
        self.genome_fasta=genome_fasta
        self.batch_size=batch_size
        self.max_jitter=max_jitter
        self.negative_sampling_ratio=negative_sampling_ratio
    
        tdb_config=get_default_config() 
        tdb_ctx=tiledb.Ctx(config=tdb_config)
        self.ctx=tdb_ctx
        self.tdb_array_name=tdb_array
        self.tdb_array=tiledb.open(tdb_array,mode='r',ctx=self.ctx)
        print("opening:"+tdb_array+" for reading...")
        print("success!")

        # using tiledb for a just a single dataset and considering counts and profile as two tasks - 
        self.num_inputs=1
        self.num_outputs=2
        self.num_tasks=2
        self.tdb_input_flank=[[input_len//2], [input_len//2]]
        self.tdb_output_flank=[[output_len//2], [output_len//2]]
        print(cts_sum_min_thresh)
        self.tdb_output_min=transform_data_type_min(["None", cts_sum_min_thresh],2)
        self.tdb_output_max=transform_data_type_max(["None", cts_sum_max_thresh],2)
        self.tdb_input_datasets=["seq"]
        self.tdb_output_datasets=[tdb_dataset, tdb_dataset]
        self.tdb_input_source_attribute=[["seq"]]
        self.tdb_output_source_attribute=[["count_bigwig_unstranded_5p"], ["count_bigwig_unstranded_5p"]]
        self.tdb_output_aggregation=[["None"], ["sum"]]
        self.tdb_output_transformation=[["None"], ["log"]]        
        self.shuffle_epoch_end=True
        self.pseudocount=1
        self.add_revcomp=add_revcomp
        self.return_coords=return_coords
            
        #  find the tdb indices that correspond to the chroms to be used 
        self.get_chrom_index_ranges()

        #get indices of input datasts to be used in training
        self.input_dataset_indices=self.get_dataset_indices(self.tdb_input_datasets)
        #get indices for outputs to be used in training
        self.output_dataset_indices=self.get_dataset_indices(self.tdb_output_datasets)

        self.map_regions_to_tdb_index()      

        print('created generator')
        
    def map_regions_to_tdb_index(self):
        self.pos_coord=[]
        self.neg_coord=[]
        self.pos_tdb_indices=[]
        self.neg_tdb_indices=[]

        if self.peak_regions is not None:
            for index,row in self.peak_regions.iterrows():
                #summit center the bed regions
                chrom=row["chr"]
                start=row["start"]
                end=row["end"]
                summit=row["summit"]
                pos=start+summit
                tdb_index=self.chrom_to_indices[chrom][0]+pos
                self.pos_coord.append([chrom,pos])
                self.pos_tdb_indices.append(tdb_index) 

        if self.nonpeak_regions is not None:
            for index,row in self.nonpeak_regions.iterrows():
                #summit center the bed regions
                chrom=row["chr"]
                start=row["start"]
                end=row["end"]
                summit=row["summit"]
                pos=start+summit
                tdb_index=self.chrom_to_indices[chrom][0]+pos
                self.neg_coord.append([chrom,pos])
                self.neg_tdb_indices.append(tdb_index) 

        subsample_neg=np.random.randint(low=0, high=len(self.neg_tdb_indices), size=(int(self.negative_sampling_ratio*len(self.neg_tdb_indices),)))
  
        self.tdb_indices = self.pos_tdb_indices+np.array(self.neg_tdb_indices)[subsample_neg].tolist()
        self.coord = self.pos_coord + np.array(self.neg_coord)[subsample_neg].tolist()

        temp = list(zip(self.coord, self.tdb_indices)) 
        random.shuffle(temp) 
        self.coord, self.tdb_indices = zip(*temp)
        
    def get_chrom_index_ranges(self):
        '''
        find tdb indices corresponding to the used chromosomes 
        '''

        num_chroms=self.tdb_array.meta['num_chroms']
        self.dataset_indices=[i for i in range(self.tdb_array.meta['num_tasks'])]
        chrom_indices=[]
        chroms=[]
        for i in range(num_chroms):
            chrom_name=self.tdb_array.meta['chrom_'+str(i)]
            chroms.append(chrom_name)
            start_index=self.tdb_array.meta['offset_'+str(i)]
            end_index=start_index+self.tdb_array.meta['size_'+str(i)]
            chrom_indices.append((start_index,end_index))   

        self.chrom_indices=chrom_indices
        self.chroms_to_use=chroms
        self.chrom_to_indices={}

        for i in range(num_chroms):
            cur_chrom=self.chroms_to_use[i]
            cur_indices=self.chrom_indices[i]
            self.chrom_to_indices[cur_chrom]=cur_indices 
    
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
    
        
    def __len__(self):
        return int(ceil(len(self.tdb_indices)/self.batch_size))
        
    def __getitem__(self,idx):

        gc.unfreeze()
        self.ref=pysam.FastaFile(self.genome_fasta)
        
        # add jitter here 
        jitter_values = np.random.randint(low=-self.max_jitter, high=self.max_jitter+1, size=(self.batch_size,))
        tdb_batch_indices=np.array(self.tdb_indices[idx*self.batch_size:(idx+1)*self.batch_size]) + jitter_values  
        
        coords=[[chrom, int(pos)+jitter_values.tolist()[i]] for i,[chrom,pos] in enumerate(self.coord[idx*self.batch_size:(idx+1)*self.batch_size])]

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
                    cur_seq=one_hot_encode(self.get_seq(coords, self.tdb_input_flank[cur_input_index][cur_input_channel_index]))
                    if cur_x is None: 
                        cur_x=cur_seq
                    else:
                        cur_x=np.concatenate((cur_x,cur_seq),axis=-1)
            
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
                cur_vals=self.get_tdb_vals(tdb_batch_indices,cur_output_channel,self.tdb_output_flank[cur_output_index][cur_output_channel_index],self.output_dataset_indices[cur_output_index][cur_output_channel_index])
                aggregate_vals=self.aggregate_vals(cur_vals,self.tdb_output_aggregation[cur_output_index][cur_output_channel_index])
                transformed_vals=self.transform_vals(aggregate_vals,self.tdb_output_transformation[cur_output_index][cur_output_channel_index])
                if cur_y is None:
                    cur_y=transformed_vals
                else:
                    cur_y=np.concatenate((cur_y,transformed_vals),axis=-1)

            #perform reverse complementation, if specified
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
        else:
            coords=None

        filtered_X,filtered_y,filtered_coords=self.remove_data_out_of_range(X,y,coords)
        
        if filtered_X[0].size==0:
            print("no data points in this batch moving to next batch")
            try:
                return self.__getitem__(idx+1)
            except:
                #we are at the last index, wrap around 
                return self.__getitem__(0)

        return (filtered_X,filtered_y,filtered_coords)
        
    def remove_data_out_of_range(self,X,y,coords=None):
        bad_indices=[]
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
        elif transformer == 'log':
            return np.log(vals+self.pseudocount)
        else:
            raise Exception("transform_vals argument must be one of None, log; you provided:"+transformer) 
    
    def aggregate_vals(self,vals,aggregator):
        if aggregator is None:
            return vals
        if aggregator == 'None':
            return vals
        elif aggregator == 'sum':
            return np.sum(vals,axis=1) 
        else:
            raise Exception("aggregate_vals argument must be one of None, sum; you provided:"+aggregator)

    
    def on_epoch_end(self):
        subsample_neg=np.random.randint(low=0, high=len(self.neg_tdb_indices), size=(int(self.negative_sampling_ratio*len(self.neg_tdb_indices),)))
        self.tdb_indices = self.pos_tdb_indices+np.array(self.neg_tdb_indices)[subsample_neg].tolist()
        self.coord = self.pos_coord + np.array(self.neg_coord)[subsample_neg].tolist()

        if self.shuffle_epoch_end==True:
            temp = list(zip(self.coord, self.tdb_indices))
            random.shuffle(temp)
            self.coord, self.tdb_indices = zip(*temp)
