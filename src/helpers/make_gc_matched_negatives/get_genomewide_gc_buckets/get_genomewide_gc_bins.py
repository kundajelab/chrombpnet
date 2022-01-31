import pandas as pd
import pysam
import argparse
import pdb 
from tqdm import tqdm

def parse_args():
    parser=argparse.ArgumentParser(description="get gc content after binning the entire genome into bins")
    parser.add_argument("-g","--genome", help="reference genome file")
    parser.add_argument("-c","--chrom_sizes", help="chromosome sizes file for reference genome (contains chr and chrom size seperated by tab)")
    parser.add_argument("-o","--output_prefix", help="output prefix path to store the gc content of binned genome")
    parser.add_argument("-f","--inputlen",type=int,default=2114, help="inputlen to use to find gc content")
    parser.add_argument("-s","--stride",type=int,default=50, help="stride to use for shifting the bins")
    return parser.parse_args()


def main():
    args=parse_args()
    ref=pysam.FastaFile(args.genome)
    chrom_sizes=pd.read_csv(args.chrom_sizes,header=None,sep='\t')
    region_dict=dict()
    for index,row in chrom_sizes.iterrows():
        chrom=row[0]
        print(chrom) 
        chrom_size=row[1]
        for bin_start in tqdm(range(0,chrom_size,args.stride)):
            bin_end=bin_start+args.inputlen
            seq=ref.fetch(chrom,bin_start,bin_end).upper()
            g=seq.count('G')
            c=seq.count('C')
            gc=g+c
            fract=round(gc/(args.inputlen),2)
            region_dict[tuple([chrom,bin_start,bin_end])]=fract

    #generate pandas df from dict
    print("making df") 
    df=pd.DataFrame.from_dict(region_dict,orient='index')
    print("made df")
    new_index=pd.MultiIndex.from_tuples(df.index, names=('CHR', 'START','END'))
    df = pd.DataFrame(df[0], new_index)
    df.to_csv(args.output_prefix+".bed",sep='\t', header=False, index=True, index_label=['CHROM','START','END'])
  
if __name__=="__main__":
    main()
