import sys 
import pandas as pd
import numpy as np 
data=pd.read_csv(sys.argv[1],header=0,sep='\t')
data['pos0']=data['position']-1
data['logratio']=np.log((data['POSTfreq']+.01)/(data['prechipfreq']+0.01))
data['logP']=-1*np.log(data['pvalue'])
bed=data[['Chr','pos0','position','logP','logratio']]
bed.to_csv(sys.argv[1]+".bed",sep='\t',header=True,index=False)


