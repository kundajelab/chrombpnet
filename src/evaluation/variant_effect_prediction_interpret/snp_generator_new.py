from tensorflow.keras.utils import Sequence
import pandas as pd
import numpy as np
import math
import pyfaidx
import one_hot

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

        cur_chrom=cur_entries.loc[0,"CHR"]
        mid = cur_entries.loc[0,"POS0"]+cur_entries.loc[1,"POS0"]
        mid = mid // 2
        print(mid, mid-cur_entries.loc[0,"POS0"], cur_entries.loc[1,"POS0"]-mid)
        left_flank_start=mid-flank_size
        left_flank_end=mid+flank_size
        sequence = str(self.genome[cur_chrom][left_flank_start:left_flank_end])

        rsid=cur_chrom+"_"+str(mid)+"_"+cur_entries.loc[0,"REF"]+"_"+cur_entries.loc[0,"ALT"]+"_"+cur_entries.loc[0,"META_DATA"]

        tt1 = mid-cur_entries.loc[0,"POS0"]
        tt2 = cur_entries.loc[1,"POS0"]-mid

        print(cur_entries.loc[0,"REF"], cur_entries.loc[1,"REF"])
        cur_ref_seq = list(sequence)
        cur_ref_seq[flank_size-tt1] = cur_entries.loc[0,"REF"]
        cur_ref_seq[flank_size+tt2] = cur_entries.loc[1,"REF"]
        cur_ref_seq=''.join(cur_ref_seq)

        cur_alt_seq = list(sequence)
        cur_alt_seq[flank_size-tt1] = cur_entries.loc[0,"ALT"]
        cur_alt_seq[flank_size+tt2] = cur_entries.loc[1,"ALT"]
        cur_alt_seq=''.join(cur_alt_seq)

        print(cur_ref_seq[flank_size-tt1-5:flank_size-tt1+1])
        print(cur_ref_seq[flank_size+tt2:flank_size+tt2+5])
        print(cur_alt_seq[flank_size-tt1-5:flank_size-tt1+1])
        print(cur_alt_seq[flank_size+tt2:flank_size+tt2+5])

        print(cur_ref_seq[flank_size-tt1], cur_ref_seq[flank_size+tt2])
        print(cur_alt_seq[flank_size-tt1], cur_alt_seq[flank_size+tt2])
           
            
        ref_seqs.append(cur_ref_seq)
        alt_seqs.append(cur_alt_seq)
        rsids.append(rsid)

        #ref_seqs.append(cur_ref_seq)
        #alt_seqs.append(cur_alt_seq)
        #rsids.append(rsid)

        return rsids, one_hot.dna_to_one_hot(ref_seqs), one_hot.dna_to_one_hot(alt_seqs)

    def __len__(self):
        return math.ceil(self.num_snps/self.batch_size)

