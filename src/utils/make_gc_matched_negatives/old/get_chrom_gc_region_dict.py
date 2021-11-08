import argparse
import pickle
from tqdm import tqdm

def parse_args():
    parser=argparse.ArgumentParser(description="generate a pickle with chrom->gc->regions")
    parser.add_argument("--input_bed",help="bed file with gc content in 4th column rounded to 2 decimals")
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.input_bed,'r')
    gc_dict={}
    index=0 
    #print(len(data))
    for line in tqdm(data):
        line=line.strip('\n') 
        index+=1
        #if index%1000000==0:
        #    print(index) 
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
    print("pickling!") 
    #pickle!
    pickle.dump( gc_dict, open( args.outf, "wb" ) )


    
if __name__=="__main__":
    main()
    
