import pandas as pd
import os
import deepdish as dd
import numpy as np

#data = pd.read_csv("model_dir_atac.csv",header=None)
#ddtpe="ATAC"
#ddtpen=ddtpe+"_PE"
#cell_types=["H1ESC" ,"HEPG2", "K562", "GM12878", "IMR90"]
#itype="counts"
##itype="profile"

#data = pd.read_csv("model_dir_dnase.csv",header=None)
data = pd.read_csv("model_dir_atac.csv",header=None)
#data = pd.read_csv("v1/model_dir_dnase_v2_interpret.csv",header=None)
#ddtpe="DNASE"
ddtpe="ATAC"
ddtpen=ddtpe+"_PE"
#ddtpen=ddtpe+"_SE"
#cell_types=["H1ESC" ,"HEPG2", "K562", "GM12878", "IMR90"]
cell_types=["GM12878", "IMR90", "H1ESC"]
#cell_types=["H1ESC"]
#cell_types=["IMR90"]
#cell_types=["HEPG2", "K562"]
#itype="counts"
itype="profile"

#h1dnase_done_jan_4_2024
odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"
for cell_type in cell_types:
	ndata = data[data[1]==cell_type].reset_index()
	#cell_type = cell_type+"_new"
	for i,r in ndata.iterrows():

		print(i,r[2])


		ppath = os.path.join(r[2],"chrombpnet_model/interpret_dnase/full_"+cell_type+"."+itype+"_scores.h5")
		print(ppath)
		if os.path.exists(ppath):
			ppath = os.path.join(r[2],"chrombpnet_model/interpret_dnase/full_"+cell_type+"."+itype+"_scores.h5")
		elif os.path.exists(os.path.join(r[2],"interpret_dnase/merged."+cell_type+"."+itype+"_scores.h5")):
			ppath = os.path.join(r[2],"interpret_dnase/merged."+cell_type+"."+itype+"_scores.h5")
		elif os.path.exists(os.path.join(odir+"/"+cell_type+"/"+r[2].split("/")[-1],"chrombpnet_model/interpret_dnase/full_"+cell_type+"."+itype+"_scores.h5")):
			ppath = os.path.join(odir+"/"+cell_type+"/"+r[2].split("/")[-1],"chrombpnet_model/interpret_dnase/full_"+cell_type+"."+itype+"_scores.h5")
			print(ppath)
		else:
			ppath = os.path.join(r[2],"chrombpnet_model/interpret_dnase/full_"+cell_type+"."+itype+"_scores.h5")
			assert(os.path.exists(ppath))
		output_prefix=ppath.replace(".h5", "_new")
		print(output_prefix)
		if os.path.isfile(output_prefix+"_compressed.h5"):
			#os.remove(output_prefix+"_compressed.h5")
			#print("removed "+output_prefix+"_compressed.h5")
			continue

		scores = dd.io.load(ppath)
		command = "python compress_deepshap.py -i "+ppath+" -o "+output_prefix
		print(command)
		os.system(command)
