#subsample regions from dl genomewide test set for use in testing SVM (using the full genomewide test set is computationally not possible)
import random
random.seed(1234) 
import argparse
import pandas as pd
import pdb 
from kerasAC.splits import * 
def parse_args():
    parser=argparse.ArgumentParser()
    parser.add_argument("--hdf_input",default="/srv/scratch/annashch/5_cell_lines_bias_correction/genomewide_labels/classificationlabels.SummitWithin200bpCenter.hdf5")
    parser.add_argument("--npos",default=50000)
    parser.add_argument("--nneg",default=50000)
    parser.add_argument("--genome",default='hg38') 
    parser.add_argument("--out_prefix",default="genomewide_svminputs/")
    parser.add_argument("--numfolds",type=int,default=10) 
    return parser.parse_args()

def main():
    args=parse_args()
    data=pd.read_hdf(args.hdf_input)
    print("loaded data") 
    colnames=data.columns
    for fold in range(args.numfolds):
        print("fold:"+str(fold))
        cur_chroms=splits[args.genome][fold]['test']

        sub_df=data.loc[cur_chroms]
        print("made sub df") 
        for colname in colnames:
            print("task:"+str(colname))
            pos=sub_df[colname]==1
            num_pos=pos.shape[0]
            pos_subset=random.sample(range(num_pos),args.npos)
            pos_index=sub_df.index[pos_subset]
            print("got pos") 
            outf=open(args.out_prefix+str(fold)+"."+colname+".positives",'w')
            for entry in pos_index:
                outf.write(entry[0]+'\t'+str(entry[1])+'\t'+str(entry[2])+'\n')
            outf.close()
            neg=sub_df[colname]==0
            num_neg=neg.shape[0]
            neg_subset=random.sample(range(num_neg),args.nneg)
            neg_index=sub_df.index[neg_subset]
            print("got neg") 
            outf=open(args.out_prefix+str(fold)+"."+colname+".negatives","w")
            for entry in neg_index:
                outf.write(entry[0]+'\t'+str(entry[1])+'\t'+str(entry[2])+'\n')
            outf.close()

            
        
    
    

if __name__=="__main__":
    main()
    
