import pandas as pd
from sklearn.metrics import average_precision_score

import sys
fold=sys.argv[1]
labels=pd.read_csv("labels."+fold+".sorted.bed",header=None,sep='\t')
preds=pd.read_csv("genomewidepredictions.IMR90."+fold+".all.sorted.bed",header=None,sep='\t')
merged=pd.concat([labels,preds],axis=1).dropna()[3]
merged.columns=['labels','preds']
cur_auprc=average_precision_score(merged['labels'],merged['preds'])
total=merged.shape[0]
pos=sum(merged['labels'])
neg=total-pos
print('IMR90\t'+fold+'\t'+str(cur_auprc)+'\t'+str(pos)+'\t'+str(neg))


