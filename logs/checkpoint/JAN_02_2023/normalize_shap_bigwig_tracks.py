import argparse
import pyBigWig
import pandas as pd
import numpy as np

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8', 'chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

chromsizes=[i.split('\t') for i in open("reference/chrom.sizes",'r').read().strip().split('\n')]
chromsizes=[tuple([i[0],int(i[1])]) for i in chromsizes]

def scale_bigwig_and_create_new(bigwig, new_pth, lower_lim, upper_lim):
        new_bigwig=pyBigWig.open(new_pth,'w')
        new_bigwig.addHeader(chromsizes)
        for val in chromsizes:
                if val[0] not in chroms:
                    continue
                print(val[0], val[1])
                old_val = np.nan_to_num(bigwig.values(val[0], 0, val[1], numpy=True))
                old_val[old_val>upper_lim] = upper_lim
                old_val[old_val<lower_lim] = lower_lim
                new_val = old_val
                new_bigwig.addEntries(val[0],0,values=new_val,span=1,step=1)
        new_bigwig.close()
        print("done")
        return

parser=argparse.ArgumentParser(description="Normalize bigwigs")
parser.add_argument("-stat","--input_stat", help="bed format - merged bed regions")
parser.add_argument("-bigwig","--input_bigwig", help="input bigwig to normalize")
parser.add_argument("-o","--output_path", help="output path")
args = parser.parse_args()

NARROWPEAK_SCHEMA = ["chr", "start", "end"] 


bw = pyBigWig.open(args.input_bigwig) 
stats = pd.read_csv(args.input_stat,sep="\t", names=["percentile", "value"])

upper_lim = float(stats[stats["percentile"]=="99.9%"]["value"])
lower_lim = float(stats[stats["percentile"]==".1%"]["value"])
print(lower_lim,upper_lim)



new_pth = args.output_path
print(new_pth)
scale_bigwig_and_create_new(bw, new_pth, lower_lim, upper_lim)
