import sys
fold=sys.argv[1]
import pandas as pd
data=pd.read_hdf("IMR90."+str(fold)+".classification.withgc.dl.pred.svmtrainset.genometestset.hdf5.labels.0")
data.to_csv("labels."+str(fold)+".bed",sep='\t',index=True,header=False)
