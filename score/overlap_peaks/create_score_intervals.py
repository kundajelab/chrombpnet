import argparse
import pandas as pd

def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--bed")
    parser.add_argument("--flank",type=int,default=1057)
    parser.add_argument("--outf")
    return parser.parse_args()

def main():
    args=parse_args()
    data=open(args.bed,'r').read().strip().split('\n')
    training_size=2*args.flank
    outf=open(args.outf,'w') 
    for region in data:
        tokens=region.split('\t')
        chrom=tokens[0]
        start_pos=int(tokens[1])
        end_pos=int(tokens[2])
        observed_region_size=end_pos-start_pos
        
        if observed_region_size < training_size:
            #can fit whole region into single training input sequence 
            center=(start_pos+end_pos)//2
            score_start=center-args.flank
            score_end=center+args.flank
            outf.write(chrom+'\t'+str(score_start)+'\t'+str(score_end)+'\n')
        else:
            #tile the region 
            score_start=start_pos
            score_end=start_pos+2*args.flank
            center=(score_start+score_end)//2
            outf.write(chrom+'\t'+str(score_start)+'\t'+str(score_end)+'\n')
            while score_end < end_pos:
                score_start=score_end+1
                score_end=score_start+2*args.flank
                center=(score_start+score_end)//2
                outf.write(chrom+'\t'+str(score_start)+'\t'+str(score_end)+'\n')
    outf.close()
    

if __name__=="__main__":
    main()
    
