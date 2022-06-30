from tensorflow.keras.utils import Sequence
import pandas as pd
import numpy as np
import math
import pyfaidx
from chrombpnet.training.utils import one_hot

class SNPGenerator(Sequence):
    def __init__(self,
                 snp_regions,
                 inputlen,
                 genome_fasta,
                 batch_size=50,
                 debug_mode_on=False):

        self.snp_regions=snp_regions
        self.num_snps=self.snp_regions.shape[0]
        self.inputlen=inputlen
        self.batch_size=batch_size
        self.genome=pyfaidx.Fasta(genome_fasta)
        self.debug_mode_on=debug_mode_on

    def __getitem__(self,idx):
        
        ref_seqs=[]
        alt_seqs=[]
        rsids=[]
            
        cur_entries=self.snp_regions.iloc[idx*self.batch_size:min([self.num_snps,(idx+1)*self.batch_size])]
        flank_size = self.inputlen // 2

        for index,entry in cur_entries.iterrows():

            cur_chrom=str(entry["CHR"])
            cur_pos=int(entry["POS0"])
            ref_snp=str(entry["REF"])
            alt_snp=str(entry["ALT"])
            meta=str(entry["META_DATA"])

            rsid=cur_chrom+"_"+str(cur_pos)+"_"+ref_snp+"_"+alt_snp+"_"+meta

            # get all regions left of snp insert locus
            left_flank_start=max([0,cur_pos-flank_size])
            left_flank_end=cur_pos
            left_flank=str(self.genome[cur_chrom][left_flank_start:left_flank_end])

            # get all regions right of snp insert locus
            right_flank_start=cur_pos+1
            right_flank_end=cur_pos+flank_size
            right_flank=str(self.genome[cur_chrom][right_flank_start:right_flank_end])

            # insert snp
            cur_ref_seq=left_flank+ref_snp+right_flank
            cur_alt_seq=left_flank+alt_snp+right_flank

            if self.debug_mode_on:
                print("CHR_POS_REF_ALT_META : " + cur_chrom+"_"+str(cur_pos)+"_"+ref_snp+"_"+alt_snp+"_"+meta + "\n")
                print("reference/alternate allele right flank : " +  right_flank + "\n")
                print("reference/alternate allele left flank : " + left_flank + "\n")
           
            if len(cur_ref_seq) != self.inputlen or len(cur_alt_seq) != self.inputlen:
                print("Exception input size is not 2114 - skipping snp")
                print("rsid (chr_pos_ref_alt): ", rsid)
                continue
            
            ref_seqs.append(cur_ref_seq)
            alt_seqs.append(cur_alt_seq)
            rsids.append(rsid)

        return rsids, one_hot.dna_to_one_hot(ref_seqs), one_hot.dna_to_one_hot(alt_seqs)

    def __len__(self):
        return math.ceil(self.num_snps/self.batch_size)

