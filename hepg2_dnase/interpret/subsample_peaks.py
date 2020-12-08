import pandas as pd
regions=pd.read_csv("HEPG2.dnase.idr.optimal_peak.narrowPeak.gz",header=None,sep='\t')
regions=regions.sample(10000,axis=0)
regions.to_csv("HEPG2.1k.dnase.idr.optimal_peak.narrowPeak",header=False,index=False,sep='\t')
