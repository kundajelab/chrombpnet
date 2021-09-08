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
                 allele_col1,
                 allele_col2,
                 flank_size,
                 ref_fasta,
                 rsid_col=None,
                 compute_gc=False,
                 batch_size=1000,
                 control_plus=None,
                 control_minus=None,
                 expand_dims=True):
        self.bed_path=bed_path
        self.bed=pd.read_csv(self.bed_path,header=0,sep='\t')
        self.num_snps=self.bed.shape[0]
        self.chrom_col=chrom_col
        self.pos_col=pos_col
        self.allele_col1=allele_col1
        self.allele_col2=allele_col2
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
        self.control_plus=control_plus
        self.control_minus=control_minus


    def __getitem__(self,idx):
        with self.lock:
            self.ref=pysam.FastaFile(self.ref_fasta)
        cur_entries=self.bed.iloc[idx*self.batch_size:min([self.num_snps,(idx+1)*self.batch_size])]
        seqs1=[]
        seqs2=[]
        gc=[]
        rsids=[]
        for index,entry in cur_entries.iterrows():
            cur_chrom=str(entry[self.chrom_col])
            if cur_chrom.startswith('chr')==False:
                cur_chrom='chr'+cur_chrom
            #print(cur_chrom, self.pos_col)
            cur_pos=entry[self.pos_col]
            left_flank_start=max([0,cur_pos-self.flank_size])
            left_flank_end=cur_pos
            snp_allele=entry[self.allele_col1]
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
            cur_seq1=left_flank+snp_allele+right_flank

            snp_allele=entry[self.allele_col2]
            try:
                if len(snp_allele)>1:
                    snp_allele=snp_allele[0]
            except:
                print(entry)
                raise Exception() 

            cur_seq2=left_flank+snp_allele+right_flank

            assert(len(cur_seq1)==2*self.flank_size)
            assert(len(cur_seq2)==2*self.flank_size)
            seqs1.append(cur_seq1)
            seqs2.append(cur_seq2)
            if self.compute_gc==True:
                cur_gc=self.compute_gc_func(cur_seq)
                gc.append(cur_gc)
            if self.rsid_col is not None:
                rsids.append(entry[self.rsid_col])
            else:
                rsids.append(index)
            index+=1 
        seqs1=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs1])
        seqs2=np.array([[ltrdict.get(x,[0,0,0,0]) for x in seq] for seq in seqs2])
        summit = cur_pos
        #print(cur_chrom, summit)
        #print(cur_chrom,summit-2500,summit+2500)
        #control_labels_plus=np.expand_dims(np.nan_to_num(self.control_plus.values(cur_chrom,summit-2500,summit+2500)),axis=1)
        #control_labels_minus=np.expand_dims(np.nan_to_num(self.control_minus.values(cur_chrom,summit-2500,summit+2500)),axis=1)
        #control_input_profile=np.expand_dims(np.concatenate((control_labels_plus,control_labels_minus),axis=-1),axis=0)

        #control_count_plus=np.expand_dims(np.expand_dims(np.array(np.log(np.sum(control_labels_plus))),axis=0),axis=1)
        #control_count_minus=np.expand_dims(np.expand_dims(np.array(np.log(np.sum(control_labels_minus))),axis=0),axis=1)    
        #control_input_count=np.concatenate((control_count_plus,control_count_minus),axis=-1)

        if self.expand_dims==True:
            seqs1=np.expand_dims(seqs1,axis=1) 
            seqs2=np.expand_dims(seqs2,axis=1) 
        if self.compute_gc==False:
            #return [rsids,seqs1,seqs2, control_input_profile, control_input_count]
            return [rsids,seqs1,seqs2]
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

