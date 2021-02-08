import pandas as pd
import sys
import pdb
task=sys.argv[1]
chrom=sys.argv[2]
data=pd.read_csv(task+'/'+task+".classificationlabels.SummitWithin200bpCenter.bed."+chrom+".gz",header=0,sep='\t')
#get all negatives
all_neg=data[data[task]==0]
all_neg.to_csv(task+'/'+task+'.'+chrom+'.universal_negatives.txt',index=False,header=True,sep='\t')

