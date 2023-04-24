import pandas as pd
import os
import deepdish as dd

data = pd.read_csv("model_dir_atac.csv",header=None)
ndata = data[data[1]=="HEPG2"].reset_index()

for i,r in ndata.iterrows():
	#print(i,r[2])
	ppath = os.path.join(r[2],"chrombpnet_model/interpret/full_HEPG2.counts_scores.h5")
	if os.path.exists(ppath):
		continue
	elif os.path.exists(os.path.join(r[2],"interpret/merged.HEPG2.counts_scores.h5")):
		ppath = os.path.join(r[2],"interpret/merged.HEPG2.counts_scores.h5")
		continue
	elif os.path.exists(os.path.join(r[2],"chrombpnet_model/interpret/full_HEPG2.counts_scores.h5")):
		continue
	else:
		print("missing interpretation", r[2], "counts")


for i,r in ndata.iterrows():
	#print(i,r[2])
	ppath = os.path.join(r[2],"chrombpnet_model/interpret/full_HEPG2.profile_scores.h5")
	if os.path.exists(ppath):
		continue
	elif os.path.exists(os.path.join(r[2],"interpret/merged.HEPG2.profile_scores.h5")):
		ppath = os.path.join(r[2],"interpret/merged.HEPG2.profile_scores.h5")
		continue
	elif os.path.exists(os.path.join(r[2],"chrombpnet_model/interpret/full_HEPG2.profile_scores.h5")):
		continue
	else:
		print("missing interpretation", r[2], "profile")
