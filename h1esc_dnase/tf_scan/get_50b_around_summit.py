import pandas as pd
import pdb 
data=pd.read_csv("H1ESC.dnase.idr.optimal_peak.narrowPeak.gz",header=None,sep='\t')
data_subset=data[[0,1,2,9]]
data_subset['summit']=data_subset[1]+data_subset[9]
data_subset['start']=data_subset['summit']-50
data_subset['end']=data_subset['summit']+50
data_subset['chr']=data_subset[0]
data_final=data_subset[['chr','start','end','summit']]
data_final.to_csv("H1ESC.idr.peaks.50bp.around.summit.bed",sep='\t',header=False,index=False) 
