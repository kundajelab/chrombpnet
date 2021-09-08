from tensorflow.keras.utils import Sequence
import pandas as pd
import numpy as np
import random
import math
import pysam
from ..util import *
import threading
import pickle
import pdb

#generate batches of SNP data with specified allele column name and flank size
class OneHotFromBedGenerator(Sequence):
    def __init__(self,
                 bed_path,
                 flank_size,
                 ref_fasta,
                 center_choice="summit",
                 summit_col_name=9,
                 batch_size=1000,
                 expand_dims=True):
        self.bed_path=bed_path
        self.bed=pd.read_csv(self.bed_path,header=0,sep='\t')
        self.num_regions=self.bed.shape[0]
        self.flank_size=flank_size
        self.batch_size=batch_size
        self.ref_fasta=ref_fasta
        self.lock=threading.Lock()
        self.expand_dims=expand_dims
        self.center_choice=center_choice
        self.summit_col_name=summit_col_name
        
    def __getitem__(self,idx):
        with self.lock:
            self.ref=pysam.FastaFile(self.ref_fasta)
        cur_entries=self.bed.iloc[idx*self.batch_size:min([self.num_regions,(idx+1)*self.batch_size])]
        seqs=[]
        coords=[]
        for index,entry in cur_entries.iterrows():
            cur_chrom=entry[0]
            if self.center_choice=="summit":
                cur_pos=entry[1]+entry[self.summit_col_name]
            else:
                assert self.center_choice=="center"
                cur_pos=entry[1]+math.floor((entry[2]-entry[1])/2)
            start_pos=max([0,cur_pos-self.flank_size])
            end_pos=start_pos+2*self.flank_size
            cur_seq=self.ref.fetch(cur_chrom,start_pos,end_pos)
            if len(cur_seq)< 2*self.flank_size:
                topad=2*self.flank_size - len(cur_seq)
                cur_seq=cur_seq+"N"*topad
            seqs.append(cur_seq)
            coords.append(tuple([cur_chrom,int(cur_pos)]))
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        if self.expand_dims==True:
            seqs=np.expand_dims(seqs,axis=1)
        print(seqs.shape)
        return coords,seqs 

    def __len__(self):
        return math.ceil(self.num_regions/self.batch_size)

