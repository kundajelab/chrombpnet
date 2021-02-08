import pickle
import pandas as pd
import sys
task=sys.argv[1]
chrom=sys.argv[2]
data=pd.read_csv(task+'/'+task+'.'+chrom+'.negative.gc.seq',header=None,sep='\t')
print("loaded")
data_dict=dict()
for index,row in data.iterrows():
    if index%1000==0:
        print(index) 
    gc=row[3]
    chrom=row[0]
    start=row[1]
    end=row[2]
    seq=row[4]
    if seq.__contains__("N"):
        continue
    if gc not in data_dict:
        data_dict[gc]=[]
    header='_'.join([str(i) for i in [chrom,start,end,gc]])
    data_dict[gc].append(header+'\n'+seq)
#pickle!
print("pickling chrom:"+str(chrom))
with open(task+'/'+task+'.'+chrom+'.negative.gc.seq.pickle','wb') as handle:
    pickle.dump(data_dict,handle)

