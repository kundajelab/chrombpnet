import pandas as pd
import pysam
import argparse
from tqdm import tqdm 

def parse_args():
    parser=argparse.ArgumentParser(description="get gc content from a bed file")
    parser.add_argument("--input_bed")
    parser.add_argument("--ref_fasta")
    parser.add_argument("--split_chroms",action="store_true",default=False)
    parser.add_argument("--out_prefix")
    parser.add_argument("--center_summit",action="store_true",default=False)
    parser.add_argument("--flank_size",type=int,default=500)
    parser.add_argument("--store_seq",action="store_true",default=False) 
    return parser.parse_args()
def get_line_narrowPeak(row,args):
    chrom=row[0]
    start=row[1]
    end=row[2] 
    if args.center_summit==True:
        summit=start+row[9]
        start=summit-args.flank_size
        end=summit+args.flank_size
    return chrom,start,end

def get_line_hdf5(index):
    chrom=index[0]
    start=index[1]
    end=index[2]
    return chrom, start,end 


def main():
    args=parse_args()
    ref=pysam.FastaFile(args.ref_fasta)
    outputs=dict()
    outf=None
    is_narrowPeak=True
    if args.input_bed.endswith('.hdf5'):
        #load as hdf5
        is_narrowPeak=False
        data=pd.read_hdf(args.input_bed,header=0,sep='\t')
    else:
        #load csv
        data=pd.read_csv(args.input_bed,header=0,sep='\t')
    print("loaded bed file")
    num_rows=str(data.shape[0])
    print("num_rows:"+num_rows) 
    cur_row=0 
    for index,row in tqdm(data.iterrows()):
        #if cur_row%1000==0:
        #    print(str(cur_row)+"/"+num_rows)
        cur_row+=1
        if is_narrowPeak is True:
            chrom,start,end=get_line_narrowPeak(row,args)
        else:
            chrom,start,end=get_line_hdf5(index)
        #extract fasta
        seq=ref.fetch(chrom,start,end).upper()
        g=seq.count('G')
        c=seq.count('C')
        gc=g+c
        gc_fract=round(gc/len(seq),2)
        if args.split_chroms is True:
            if chrom not in outputs:
                outputs[chrom]=open(args.out_prefix+'.'+chrom,'w')
                print("created:"+str(args.out_prefix+'.'+chrom))
            outputs[chrom].write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(gc_fract))
            if args.store_seq is True:
                outputs[chrom].write('\t'+seq+'\n')
            else:
                outputs[chrom].write('\n')
        else:
            if outf is None:
                outf=open(args.out_prefix,'w')
                print("created:"+str(args.out_prefix))
            outf.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(gc_fract))
            if args.store_seq is True:
                outf.write('\t'+seq+'\n')
            else:
                outf.write('\n')
    #close files
    if args.split_chroms is True:
        for chrom in outputs:
            outputs[chrom].close()
    else:
        outf.close()
        
if __name__=="__main__":
    main()
