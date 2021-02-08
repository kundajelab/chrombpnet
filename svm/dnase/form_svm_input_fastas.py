import random
random.seed(1234)
import pickle
import argparse
import pysam
import pandas as pd


def parse_args():
    parser=argparse.ArgumentParser(description="form svm inputs")
    parser.add_argument("--outf",nargs="+")
    parser.add_argument("--overwrite_outf",action="store_true",default=False)
    parser.add_argument("--chrom")
    parser.add_argument("--neg_pickle")
    parser.add_argument("--ref_fasta")
    parser.add_argument("--peaks",nargs="+") 
    return parser.parse_args()

def scale_gc(cur_gc):
    if random.random()>0.5:
        cur_gc+=0.01
    else:
        cur_gc-=0.01
    cur_gc=round(cur_gc,2)
    if cur_gc <=0:
        cur_gc+=0.01
    if cur_gc>=1:
        cur_gc-=0.01
    assert cur_gc >=0
    assert cur_gc <=1
    return cur_gc 

def adjust_gc(cur_gc,negatives,used_negatives):
    #verify that cur_gc is in negatives dict
    if cur_gc not in used_negatives:
        used_negatives[cur_gc]=dict()
    while (cur_gc not in negatives) or (len(used_negatives[cur_gc])>=len(negatives[cur_gc])):
        #all options for this gc value have been used up. Adjust the gc
        cur_gc=scale_gc(cur_gc)
        if cur_gc not in used_negatives:
            used_negatives[cur_gc]=dict()
    return cur_gc,used_negatives 

        

def main():
    args=parse_args()
    ref=pysam.FastaFile(args.ref_fasta)
    negatives=pickle.load(open(args.neg_pickle, "rb"))
    for cur_key in negatives:
        rounded=round(cur_key,2)
        negatives[rounded]=negatives[cur_key]
    used_negatives=dict()
    for i in range(len(args.peaks)):
        cur_peaks=args.peaks[i]
        cur_outf=args.outf[i]
        print("cur_peaks:"+cur_peaks)
        print("cur_outf:"+cur_outf) 
        if args.overwrite_outf is True:
            outf_pos=open(cur_outf+'.positives','w')
            outf_neg=open(cur_outf+'.negatives','w')
        else:
            #open in append mode for the current chromosome 
            outf_pos=open(cur_outf+'.positives','a')
            outf_neg=open(cur_outf+'.negatives','a')
        #chrom->start->end->gc->seq
        cur_peaks=pd.read_csv(cur_peaks,header=None,sep='\t')
        for index,row in cur_peaks.iterrows():
            chrom=row[0]
            if chrom!=args.chrom:
                continue
            cur_gc,used_negatives=adjust_gc(row[3],negatives,used_negatives)
            start=row[1]
            end=row[2]
            seq=row[4]
            header='_'.join([str(i) for i in [chrom,start,end,row[3]]])
            
            #get matched negative sequence
            num_candidates=len(negatives[cur_gc])
            rand_neg_index=random.randint(0,num_candidates-1)
            while rand_neg_index in used_negatives[cur_gc]:
                cur_gc,used_negatives=adjust_gc(cur_gc,negatives,used_negatives)
                num_candidates=len(negatives[cur_gc])
                rand_neg_index=random.randint(0,num_candidates-1)
            #make sure we don't sample this value again
            used_negatives[cur_gc][rand_neg_index]=1
            cur_neg=negatives[cur_gc][rand_neg_index]            
            outf_pos.write('>'+header+'\n'+seq+'\n')
            outf_neg.write('>'+cur_neg+'\n')
        outf_pos.close()
        outf_neg.close() 

if __name__=="__main__":
    main()
    
    
