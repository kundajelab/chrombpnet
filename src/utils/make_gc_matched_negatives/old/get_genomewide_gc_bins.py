import pandas as pd
import pysam
import argparse
import pdb 
from tqdm import tqdm

def parse_args():
    parser=argparse.ArgumentParser(description="get gc content from a bed file")
    parser.add_argument("--chrom_sizes")
    parser.add_argument("--ref_fasta")
    parser.add_argument("--out_prefix")
    parser.add_argument("--region_size",type=int,default=1000)
    parser.add_argument("--stride",type=int,default=50)
    parser.add_argument("--output_format",choices=['tsv','hdf5'],help="store output track as either a .tsv or an .hdf5 file")
    return parser.parse_args()


def main():
    args=parse_args()
    ref=pysam.FastaFile(args.ref_fasta)
    chrom_sizes=pd.read_csv(args.chrom_sizes,header=None,sep='\t')
    region_dict=dict()
    for index,row in chrom_sizes.iterrows():
        chrom=row[0]
        print(chrom) 
        chrom_size=row[1]
        for bin_start in tqdm(range(0,chrom_size,args.stride)):
            #if bin_start%1000000==0:
            #    print(str(bin_start))
            bin_end=bin_start+2*args.region_size
            seq=ref.fetch(chrom,bin_start,bin_end).upper()
            g=seq.count('G')
            c=seq.count('C')
            gc=g+c
            fract=round(gc/(2*args.region_size),2)
            region_dict[tuple([chrom,bin_start,bin_end])]=fract
    #generate pandas df from dict
    print("making df") 
    df=pd.DataFrame.from_dict(region_dict,orient='index')
    print("made df")
    new_index=pd.MultiIndex.from_tuples(df.index, names=('CHR', 'START','END'))
    df = pd.DataFrame(df[0], new_index)
    if args.output_format=="tsv":
        df.to_csv(args.out_prefix+".tsv",sep='\t', header=True, index=True, index_label=['CHROM','START','END'])
    else:
        assert args.output_format=="hdf5"
        df.to_hdf(args.out_prefix+".hdf5",key='data',mode='w',append=False,format='table',min_itemsize=30)                
        
if __name__=="__main__":
    main()
