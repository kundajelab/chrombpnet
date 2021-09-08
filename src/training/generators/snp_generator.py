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
class SNPGenerator(Sequence):
    def __init__(self,
                 bed_path,
                 chrom_col,
                 pos_col,
                 allele_col,
                 flank_size,
                 ref_fasta,
                 rsid_col=None,
                 compute_gc=False,
                 batch_size=1000,
                 expand_dims=True):
        self.bed_path=bed_path
        self.bed=pd.read_csv(self.bed_path,header=0,sep='\t')
        self.num_snps=self.bed.shape[0]
        self.chrom_col=chrom_col
        self.pos_col=pos_col
        self.allele_col=allele_col
        self.flank_size=flank_size
        self.rsid_col=rsid_col
        self.compute_gc=compute_gc
        self.batch_size=batch_size
        self.ref_fasta=ref_fasta
        self.batch_size=batch_size
        self.ref_fasta=ref_fasta
        self.batch_size=batch_size
        self.lock=threading.Lock()
        self.expand_dims=expand_dims


    def __getitem__(self,idx):
        with self.lock:
            self.ref=pysam.FastaFile(self.ref_fasta)
        cur_entries=self.bed.iloc[idx*self.batch_size:min([self.num_snps,(idx+1)*self.batch_size])]
        seqs=[]
        gc=[]
        rsids=[]
        for index,entry in cur_entries.iterrows():
            cur_chrom=str(entry[self.chrom_col])
            if cur_chrom.startswith('chr')==False:
                cur_chrom='chr'+cur_chrom
            cur_pos=entry[self.pos_col]
            left_flank_start=max([0,cur_pos-self.flank_size])
            left_flank_end=cur_pos
            snp_allele=entry[self.allele_col]
            try:
                if len(snp_allele)>1:
                    snp_allele=snp_allele[0]
            except:
                print(entry)
                raise Exception() 
            right_flank_start=cur_pos+1
            right_flank_end=cur_pos+self.flank_size
            left_flank=self.ref.fetch(cur_chrom,left_flank_start,left_flank_end)
            right_flank=self.ref.fetch(cur_chrom,right_flank_start,right_flank_end)
            cur_seq=left_flank+snp_allele+right_flank
            #print(snp_allele, entry[self.rsid_col])
            if len(cur_seq) != 2*self.flank_size:
                print(cur_chrom,cur_pos)
                print(self.ref.fetch(cur_chrom,left_flank_start,left_flank_end+2))
                print(entry[self.rsid_col])
                print(cur_seq)
                print(len(cur_seq))
            
            assert(len(cur_seq)==2*self.flank_size)
            seqs.append(cur_seq)
            if self.compute_gc==True:
                cur_gc=self.compute_gc_func(cur_seq)
                gc.append(cur_gc)
            if self.rsid_col is not None:
                rsids.append(entry[self.rsid_col])
            else:
                rsids.append(index)
            index+=1 
        seqs=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs])
        if self.expand_dims==True:
            seqs=np.expand_dims(seqs,axis=1) 
        if self.compute_gc==False:
            return [rsids,seqs]
        else:
            gc=np.expand_dims(np.asarray(gc),axis=1)
            return [rsids,[seqs,gc]]

    def compute_gc_func(self,seq):
        seq=seq.lower()
        g_count=seq.count('g')
        c_count=seq.count('c')
        seq_len=len(seq)
        return (g_count+c_count)/seq_len 

    def __len__(self):
        return math.ceil(self.num_snps/self.batch_size)

