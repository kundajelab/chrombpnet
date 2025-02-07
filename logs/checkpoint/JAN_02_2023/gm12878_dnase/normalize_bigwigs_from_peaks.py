import argparse
import pyBigWig
import pandas as pd
import numpy as np

chroms = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8', 'chr9','chr10', 'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

chromsizes=[i.split('\t') for i in open("/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes",'r').read().strip().split('\n')]
chromsizes=[tuple([i[0],int(i[1])]) for i in chromsizes]

def scale_bigwig_and_create_new(bigwig, new_pth, scale):
        new_bigwig=pyBigWig.open(new_pth,'w')
        new_bigwig.addHeader(chromsizes)
        print(scale)
        for val in chromsizes:
                if val[0] not in chroms:
                    continue
                print(val[0], val[1])
                old_val = np.nan_to_num(bigwig.values(val[0], 0, val[1], numpy=True))
                new_val = old_val/scale
                new_bigwig.addEntries(val[0],0,values=new_val,span=1,step=1)
        new_bigwig.close()
        print("done")
        return

parser=argparse.ArgumentParser(description="Normalize bigwigs")
parser.add_argument("-bed","--input_bed", help="bed format - merged bed regions")
parser.add_argument("-bigwig","--input_bigwig", help="input bigwig to normalize")
parser.add_argument("-o","--output_path", help="output path")
args = parser.parse_args()

NARROWPEAK_SCHEMA = ["chr", "start", "end"] 


bw = pyBigWig.open(args.input_bigwig) 
peaks = pd.read_csv(args.input_bed,sep="\t", names=NARROWPEAK_SCHEMA)
print(peaks.head())

signal_in_peaks = []
num_bases = []
for i,r in peaks.iterrows():
	vals = np.nan_to_num(bw.values(r["chr"], r['start'], r['end']))
	signal_in_peaks.append(np.sum(vals))
	num_bases.append(vals.shape[0])
	#print(vals.shape[0])

scale = np.sum(signal_in_peaks)/np.sum(num_bases)
print(scale)
print(np.sum(signal_in_peaks))
print(np.sum(num_bases))


new_pth = args.output_path
f = open(new_pth, "w")
f.write(str(scale))
f.close()

#print(new_pth)

#expected_signal = np.mean(mean_signal_in_peaks)
#scale_bigwig_and_create_new(bw, new_pth, expected_signal)
#scale_bigwig_and_create_new(bw, new_pth, scale)
