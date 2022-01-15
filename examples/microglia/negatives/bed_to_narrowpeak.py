import sys
import pandas as pd
data=pd.read_csv(sys.argv[1],header=None,sep='\t')
data[3]='.'
data[4]='.'
data[5]='.'
data[6]='.'
data[7]='.'
data[8]='.'
data[9]=round(0.5*(data[2]-data[1]))
data[9]=data[9].astype(int)
data.to_csv(sys.argv[1]+".narrowPeak",header=False,index=False,sep='\t')

