#converte an hdf5 data frame test set for neural network into svm input fasta file
import sys
import pandas as pd
import argparse
import pdb 
def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--hdf5")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    num_inputs=len(args.hdf5)
    num_outputs=len(args.outf)
    #assert num_inputs==num_outputs
    outf=open(args.outf,'w')
    print("reading hdf5") 
    cur_hdf5=pd.read_hdf(args.hdf5).index
    print("loaded hdf5")
    i=0
    for entry in cur_hdf5:
        i+=1
        if i%1000==0:
            print(i)
        chrom=entry[0]
        startpos=entry[1]
        endpos=entry[2]
        outf.write(chrom+'\t'+str(startpos)+'\t'+str(endpos)+'\n')
        

if __name__=="__main__":
    main()
    
