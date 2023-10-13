import pandas as pd
import os
import deepdish as dd
import pyBigWig as pw
import numpy as np

data = pd.read_csv("model_dir_atac.csv",header=None)
ndata = data[data[1]=="K562"].reset_index()

#merged.K562.counts.bw
#K562_w_bias.bw
#full_K562_w_bias.bw

file_inti="full_K562_w_bias.bw"
file_int="K562_w_bias.bw"

output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/merge_folds/" + file_inti

fname="/mnt/lab_data2/anusri/chrombpnet/reference/chrom.sizes"
with open(fname) as f:
	gs = [x.strip().split('\t') for x in f]
gs = [(x[0], int(x[1])) for x in gs if len(x)==2]

bw_out = pw.open(output_prefix, 'w')
bw_out.addHeader(gs)


bigwigs = []
for i,r in ndata.iterrows():
	print(i,r[2])

	ppath = os.path.join(r[2],"chrombpnet_model/interpret/"+file_int)
	if os.path.exists(ppath):
		bw = pw.open(ppath)
	elif os.path.exists(os.path.join(r[2],"interpret/"+file_int)):
		ppath = os.path.join(r[2],"interpret/"+file_int)
		bw = pw.open(ppath)

	else:
		ppath = os.path.join(r[2],"chrombpnet_model/interpret/"+file_inti)
		bw = pw.open(ppath)
        
	bigwigs.append(bw)

for val in gs:
	avg_bw = None
	for bw in bigwigs:
		sizei = int(val[1]) 
		print(val[0],sizei)
		bw_val = np.nan_to_num(bw.values(val[0],0,sizei))
		assert(bw_val.shape[0]==sizei)
		numz = sum(bw_val==0)

		if avg_bw is None:
			avg_bw=bw_val
			numzi = numz
		else:
			avg_bw+=bw_val
			print(numzi,numz)
			assert(abs(numzi-numz)<1000)

	avg_bw=avg_bw/5
	bw_out.addEntries([val[0]]*sizei,list(range(0,sizei)), ends=list(range(1,sizei+1)), values=avg_bw)

bw.close()

