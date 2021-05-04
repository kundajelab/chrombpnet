from os import listdir
from os.path import isfile, join
import pandas as pd
import pdb 
peak_dir='/oak/stanford/groups/akundaje/projects/chrombpnet/bigwigs_unplugged_bias/idr'
top_n=100000
peak_files=[f for f in listdir(peak_dir) if join(peak_dir,f).endswith('narrowPeak')]
for peak_file in peak_files:
    data=pd.read_csv(join(peak_dir,peak_file),header=None,sep='\t')
    nrows=data.shape[0]
    last_row=min([nrows,top_n])
    top_peaks=data.sort_values(by=[8],ascending=False)[1:last_row]
    top_peaks.to_csv('top100.'+peak_file,index=False,header=False,sep='\t')
    
    
