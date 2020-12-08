import pandas as pd
regions=pd.read_csv("H1ESC.dnase.idr.optimal_peak.narrowPeak.gz",header=None,sep='\t')
regions=regions.sample(10000,axis=0)
regions.to_csv("H1ESC.1k.dnase.idr.optimal_peak.narrowPeak",header=False,index=False,sep='\t')
