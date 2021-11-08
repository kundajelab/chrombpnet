import pysam
import argparse
from tqdm import tqdm 
import pandas as pd

def parse_args():
    parser=argparse.ArgumentParser(description="get gc content from a foreground bed file")
    parser.add_argument("-i","--input_bed", help="bed file in narrow peak format - we will find gc content of these regions centered on the summit")
    parser.add_argument("-g", "--ref_fasta", help="reference genome fasta")
    parser.add_argument("-o", "--out_prefix", help="output file prefix for storing gc-content values of given foreground bed")
    parser.add_argument("-f","--flank_size",type=int,default=1057, help="flank size to use to fine gc-content")
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

        # calculate gc content when centered at summit
        seq=ref.fetch(chrom,start,end).upper()
        g=seq.count('G')
        c=seq.count('C')
        gc=g+c
        gc_fract=round(gc/len(seq),2)                
        outf.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(gc_fract)+"\n")
    outf.close()
        
if __name__=="__main__":
    main()
