import pandas as pd
import random
import argparse
import pysam

def parse_args():
    parser=argparse.ArgumentParser(description='select 10k regions for gkmexplain on canonical cell types')
    parser.add_argument("-narrowPeak")
    parser.add_argument("-ref_fasta") 
    parser.add_argument("-n_to_sample",type=int)
    parser.add_argument("-outf")
    parser.add_argument("-flank",type=int)
    return parser.parse_args()

    
def main():
    args=parse_args()
    ref=pysam.FastaFile(args.ref_fasta)
    peaks=pd.read_csv(args.narrowPeak,header=None,sep='\t')
    #subsample to specific number of peaks
    subset=peaks.sample(n=args.n_to_sample,axis=0,random_state=1234)
    outf=open(args.outf,'w') 
    for index,row in subset.iterrows():
        chrom=row[0]
        start_pos=row[1]
        summit=row[9]
        summit=start_pos+summit
        first=summit - args.flank 
        last=summit + args.flank 
        seq_name=">"+str(chrom)+":"+str(summit)
        seq=ref.fetch(chrom,first,last)
        assert len(seq)==2*args.flank
        outf.write(seq_name+'\n'+seq+'\n')
    outf.write('\n')
    outf.close()
        
if __name__=="__main__":
    main()
    

    
