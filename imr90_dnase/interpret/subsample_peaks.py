import pandas as pd
regions=pd.read_csv("IMR90.dnase.idr.optimal_peak.narrowPeak.gz",header=None,sep='\t')
regions=regions.sample(10000,axis=0)
regions.to_csv("IMR90.1k.dnase.idr.optimal_peak.narrowPeak",header=False,index=False,sep='\t')
