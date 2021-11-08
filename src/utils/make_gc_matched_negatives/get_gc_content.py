import pandas as pd
import pysam
import argparse
from tqdm import tqdm 

def parse_args():
    parser=argparse.ArgumentParser(description="get gc content from a foreground bed file")
    parser.add_argument("--input_bed")
    parser.add_argument("--ref_fasta")
    parser.add_argument("--out_prefix")
    parser.add_argument("--flank_size",type=int,default=500)
    return parser.parse_args()

def main():
    args=parse_args()
    ref=pysam.FastaFile(args.ref_fasta)
    data=pd.read_csv(args.input_bed,header=0,sep='\t')

    num_rows=str(data.shape[0])
    print("num_rows:"+num_rows) 

    outf=open(args.out_prefix,'w')

    for index,row in tqdm(data.iterrows()):
        chrom=row[0]
        start=row[1]
        end=row[2] 

        summit=start+row[9]
        start=summit-args.flank_size
        end=summit+args.flank_size

        seq=ref.fetch(chrom,start,end).upper()
        g=seq.count('G')
        c=seq.count('C')
        gc=g+c
        gc_fract=round(gc/len(seq),2)                
        outf.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(gc_fract)+"\n")
    outf.close()
        
if __name__=="__main__":
    main()
