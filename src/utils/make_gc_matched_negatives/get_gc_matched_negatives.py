import argparse
import pickle
from tqdm import tqdm

def parse_args():
    parser=argparse.ArgumentParser(description="generate a bed file for negatives which are gc-matched with foreground")
    parser.add_argument("--candidate_negatives",help="bed file with gc content in 4th column rounded to 2 decimals")
    parser.add_argument("--foreground_gc_bed")
    parser.add_argument("--out_prefix")
    return parser.parse_args()

def make_gc_dict():
    """
    Imports the TF-MoDISco hits as a single Pandas DataFrame.
    The `key` column is the name of the originating PFM, and `peak_index` is the
    index of the peak file from which it was originally found.
    """
    data=open(args.candidate_negatives,'r')
    gc_dict={}
    index=0 
    for line in tqdm(data):
        line=line.strip('\n') 
        index+=1
        tokens=line.split('\t')
        chrom=tokens[0]
        gc=float(tokens[-1])
        start=tokens[1]
        end=tokens[2]
        if chrom not in gc_dict:
            gc_dict[chrom]={}
        if gc not in gc_dict[chrom]:
            gc_dict[chrom][gc]=[(chrom,start,end)]
        else:
            gc_dict[chrom][gc].append([(chrom,start,end)])
        return gc_dict

def scale_gc(cur_gc):
    if random.random()>0.5:
        cur_gc+=0.01
    else:
        cur_gc-=0.01
    cur_gc=round(cur_gc,2)
    if cur_gc<=0:
        cur_gc+=0.01
    if cur_gc>=1:
        cur_gc-=0.01
    assert cur_gc >=0
    assert cur_gc <=1
    return cur_gc 

def adjust_gc(chrom,cur_gc,negatives,used_negatives):
    #verify that cur_gc is in negatives dict
    if chrom  not in used_negatives:
        used_negatives[chrom]={}

    if cur_gc not in used_negatives[chrom]:
        used_negatives[chrom][cur_gc]=[]

    while (cur_gc not in negatives[chrom]) or (len(used_negatives[chrom][cur_gc])>=len(negatives[chrom][cur_gc])):
        #all options for this gc value have been used up. Adjust the gc
        cur_gc=scale_gc(cur_gc)
        if cur_gc not in used_negatives[chrom]:
            used_negatives[chrom][cur_gc]=dict()
    return cur_gc,used_negatives 

        
def main():

    args=parse_args()
    negatives=make_gc_dict(args)
    used_negatives=dict()
    cur_peaks=pd.read_csv(cur_peaks,header=None,sep='\t')
 
    for index,row in tqdm(cur_peaks.iterrows()): 
        chrom=row[0]
        start=row[1]
        end=row[2]
        gc_value=row[3]

        cur_gc,used_negatives=adjust_gc(chrom,gc_value,negatives,used_negatives)
        
        header='_'.join([str(i) for i in [chrom,start,end,row[3]]])
        #get matched negative sequence

        num_candidates=len(negatives[chrom][cur_gc])
        rand_neg_index=random.randint(0,num_candidates-1)
        while rand_neg_index in used_negatives[chrom][cur_gc]:
            cur_gc,used_negatives=adjust_gc(chrom,cur_gc,negatives,used_negatives)
            num_candidates=len(negatives[chrom][cur_gc])
            rand_neg_index=random.randint(0,num_candidates-1)

        used_negatives[chrom][cur_gc].append(rand_neg_index)
        neg_tuple=negatives[chrom][cur_gc][rand_neg_index][0]
        neg_chrom=neg_tuple[0]
        neg_start=neg_tuple[1]
        neg_end=neg_tuple[2]
        
        neg_header='_'.join([str(i) for i in [neg_chrom,neg_start,neg_end,cur_gc]])

def main():
    args=parse_args()

    
if __name__=="__main__":
    main()
    
