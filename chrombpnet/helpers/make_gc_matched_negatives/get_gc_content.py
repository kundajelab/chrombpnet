import pyfaidx
import argparse
from tqdm import tqdm 
import pandas as pd

def parse_args():
    parser=argparse.ArgumentParser(description="get gc content from a foreground bed file")
    parser.add_argument("-i","--input_bed", help="bed file in narrow peak format - we will find gc content of these regions centered on the summit")
    parser.add_argument("-c","--chrom_sizes",type=str, required=True, help="TSV file with chromosome name in first column and size in the second column")
    parser.add_argument("-g", "--genome", help="reference genome fasta")
    parser.add_argument("-op", "--output_prefix", help="output file prefix for storing gc-content values of given foreground bed")
    parser.add_argument("-il","--inputlen",type=int,default=2114, help="inputlen to use to find gc-content")
    return parser.parse_args()

def main(args):
    chrom_sizes_dict = {line.strip().split("\t")[0]:int(line.strip().split("\t")[1]) for line in open(args.chrom_sizes).readlines()}
    ref=pyfaidx.Fasta(args.genome)
    data=pd.read_csv(args.input_bed,header=None,sep='\t')
    assert(args.inputlen%2 == 0) # for symmtery

    num_rows=str(data.shape[0])
    print("num_rows:"+num_rows) 

    outf=open(args.output_prefix+".bed",'w')
    filtered_points=0
    for index,row in tqdm(data.iterrows()):
        chrom=row[0]
        start=row[1]
        end=row[2] 

        summit=start+row[9]
        start=summit-args.inputlen//2
        end=summit+args.inputlen//2

        if start < 0:
            filtered_points+=1
            continue
        if end > chrom_sizes_dict[chrom]:
            filtered_points+=1
            continue

        # calculate gc content when centered at summit
        seq=str(ref[chrom][start:end]).upper()
        g=seq.count('G')
        c=seq.count('C')
        gc=g+c
        gc_fract=round(gc/len(seq),2)                
        outf.write(chrom+'\t'+str(start)+'\t'+str(end)+'\t'+str(gc_fract)+"\n")
    outf.close()

    print("Number of regions filtered because inputlen sequence cannot be constructed: " + str(filtered_points))
    print("Percentage of regions filtered " + str(round(filtered_points*100.0/data.shape[0],3)) + "%" )
    if round(filtered_points*100.0/data.shape[0],3) > 25:
    	print("WARNING: If percentage of regions filtered is high (>25%) - your genome is very small - consider using a reduced input/output length for your genome")
        
if __name__=="__main__":
    args=parse_args()
    main(args)
