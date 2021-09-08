#predictions from tiledb
from tensorflow.keras.utils import Sequence
import pandas as pd
import numpy as np
import random
import threading 
from random import shuffle 
import math
from math import ceil, floor
import pysam
from ..util import *
from .tiledb_generator import * 
import tiledb
import pdb

class TiledbPredictGenerator(TiledbGenerator):
    def __init__(self,
                 ref_fasta,
                 batch_size,
                 tdb_array,
                 tdb_partition_attribute_for_upsample,
                 tdb_partition_thresh_for_upsample,
                 tdb_partition_datasets_for_upsample,
                 upsample_ratio,
                 tdb_input_source_attribute,
                 tdb_input_flank,
                 tdb_output_source_attribute,
                 tdb_output_flank,
                 tdb_ambig_attribute,
                 num_inputs,
                 num_outputs,
                 tdb_input_min=None,
                 tdb_input_max=None,
                 tdb_output_min=None,
                 tdb_output_max=None,
                 tdb_input_aggregation=None,
                 tdb_input_transformation=None,
                 tdb_output_aggregation=None,
                 tdb_output_transformation=None,
                 chroms=None,
                 chrom_sizes=None,
                 tdb_input_datasets=None,
                 tdb_output_datasets=None,
                 pseudocount=0.001,
                 tiledb_stride=1,
                 tdb_config=None,
                 tdb_ctx=None,
                 num_threads=1,
                 bed_regions=None,
                 bed_regions_center=None,
                 add_revcomp=False):
        
        TiledbGenerator.__init__(self,          
                                 ref_fasta=ref_fasta,
                                 batch_size=batch_size,
                                 tdb_array=tdb_array,
                                 tdb_partition_attribute_for_upsample=tdb_partition_attribute_for_upsample,
                                 tdb_partition_thresh_for_upsample=tdb_partition_thresh_for_upsample,
                                 tdb_partition_datasets_for_upsample=tdb_partition_datasets_for_upsample,
                                 upsample_ratio=upsample_ratio,
                                 tdb_input_source_attribute=tdb_input_source_attribute,
                                 tdb_input_flank=tdb_input_flank,
                                 tdb_input_min=tdb_input_min,
                                 tdb_input_max=tdb_input_max,
                                 tdb_output_min=tdb_output_min,
                                 tdb_output_max=tdb_output_max,
                                 tdb_output_source_attribute=tdb_output_source_attribute,
                                 tdb_output_flank=tdb_output_flank,
                                 num_inputs=num_inputs,
                                 num_outputs=num_outputs,
                                 tdb_input_aggregation=tdb_input_aggregation,
                                 tdb_input_transformation=tdb_input_transformation,
                                 tdb_output_aggregation=tdb_output_aggregation,
                                 tdb_output_transformation=tdb_output_transformation,
                                 tdb_ambig_attribute=tdb_ambig_attribute,
                                 chroms=chroms,
                                 chrom_sizes=chrom_sizes,
                                 shuffle_epoch_start=False,
                                 shuffle_epoch_end=False,
                                 pseudocount=pseudocount,
                                 add_revcomp=add_revcomp,
                                 return_coords=True,
                                 tdb_config=tdb_config,
                                 tdb_ctx=tdb_ctx,
                                 tdb_input_datasets=tdb_input_datasets,
                                 tdb_output_datasets=tdb_output_datasets,
                                 num_threads=num_threads,
                                 bed_regions=bed_regions,
                                 bed_regions_center=bed_regions_center)
        self.tiledb_stride=tiledb_stride
        self.idx_to_tdb_index=None
        print("created predict generator")
        
    def precompute_idx_to_tdb_index(self):
        self.idx_to_tdb_index={}
        cur_chrom_index=0
        #print(self.chrom_indices[cur_chrom_index])
        #print(self.tdb_input_flank)
        #print(self.chrom_indices[cur_chrom_index][0],self.tdb_input_flank[0])
        next_coord=self.chrom_indices[cur_chrom_index][0]+self.tdb_input_flank[0][0] #first tdb coord for chromosome,this is next batch's start coord 
        cur_chrom_end=self.chrom_indices[cur_chrom_index][1]-self.tdb_input_flank[0][0] #final tdb coordinate for this chromsome
        print("mapping idx to tiledb indices") 
        for idx in range(0,self.__len__()):
            self.idx_to_tdb_index[idx]=[]
            for i in range(self.batch_size): 
                if next_coord > cur_chrom_end:
                    try:
                        cur_chrom_index+=1
                        next_coord=self.chrom_indices[cur_chrom_index][0]+self.tdb_input_flank[0][0] 
                        cur_chrom_end=self.chrom_indices[cur_chrom_index][1]-self.tdb_input_flank[0][0] 
                    except:
                        #we are in the last batch, which may not be full-sized, we cannot increment to the next chromosome. 
                        continue 
                self.idx_to_tdb_index[idx].append(next_coord)
                next_coord+=self.tiledb_stride
        print("mapping of idx to tiledb indices completed")
        print("MAX idx:"+str(idx))
    def get_tdb_indices_for_batch(self,idx):
        if self.bed_regions is not None:
            return self.tdb_indices[idx*self.batch_size:(idx+1)*self.batch_size]
        elif len(self.upsampled_indices)>0:
            #use the upsampled indices 
            upsampled_batch_start=idx*self.upsampled_batch_size
            upsampled_batch_end=min([upsampled_batch_start+self.upsampled_batch_size,self.upsampled_indices.shape[0]])
            upsampled_batch_indices=self.upsampled_indices[upsampled_batch_start:upsampled_batch_end]
            return upsampled_batch_indices
        else:
            #not upsampling, going through the test chromosomes with specified stride value
            #generate mapping of idx values to tdb indices 
            if self.idx_to_tdb_index is None:
                self.precompute_idx_to_tdb_index()
            batch_indices=self.idx_to_tdb_index[idx]
            return batch_indices
    
    def __len__(self):
        #we have an explict set of regions
        if self.bed_regions is not None:
            return int(ceil(len(self.tdb_indices)/self.batch_size))
        elif len(self.upsampled_indices) is 0: 
            return int(ceil(self.num_indices/(self.batch_size*self.tiledb_stride)))
        else:
            return int(ceil(self.upsampled_indices.shape[0] /self.batch_size))

    def on_epoch_end(self):
        pass

