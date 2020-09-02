import tiledb
import numpy as np 
import pdb 
datapath="/mnt/lab_data2/projects/adpd/adpd_tiledb_pseudobulk/Cluster24.chr1"
data=tiledb.DenseArray(datapath,mode='r')
fc_bigwig=data[:]['fc_bigwig']
idr_peak=data[:]['idr_peak']
overlap_peak=data[:]['overlap_peak']
pval_bigwig=data[:]['pval_bigwig']
ambig_peak=data[:]['ambig_peak']
#ambig_peak=data[:]['count_bigwig_plus_5p']
#ambig_peak=data[:]['count_bigwig_minus_5p']
print(np.nanmax(overlap_peak))
print(np.nanmax(idr_peak))
print(np.nanmax(ambig_peak))
print(np.nanmax(pval_bigwig))
print(np.nanmax(fc_bigwig))

pdb.set_trace() 
