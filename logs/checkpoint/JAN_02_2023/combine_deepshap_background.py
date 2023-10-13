import pandas as pd
import os
import deepdish as dd
import numpy as np
from tqdm import tqdm

data = pd.read_csv("model_dir_atac.csv",header=None)
ddtpe="ATAC"
ddtpen=ddtpe+"_PE"
#cell_types=["HEPG2", "K562", "GM12878", "H1ESC", "IMR90"]
#cell_types=["K562", "GM12878", "H1ESC", "IMR90", "HEPG2"]
cell_types=["IMR90", "HEPG2"]
#cell_types=["HEPG2"]
itype="profile"


#data = pd.read_csv("model_dir_dnase.csv",header=None)


NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

				
	
for cell_type in cell_types:
	ndata = data[data[1]==cell_type].reset_index()
	bed_of_interest = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"+ddtpen+"/"+cell_type+"/negatives_data/negatives_subsample.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
	one_hots=None
	for i,r in ndata.iterrows():
		print(i,r[2])

		beds_path = os.path.join(r[2],"background_interpret/"+cell_type+".interpreted_regions_v2.bed")
		if os.path.exists(beds_path):
			beds = pd.read_csv(beds_path, sep="\t", header=None)

		print(beds.head())
		beds["key"] = beds[0] + "_" + beds[1].astype(str) + "_" + beds[2].astype(str) + "_" + + beds[9].astype(str)

		ppath = os.path.join(r[2],"background_interpret/"+cell_type+".profile_scores.h5")
		if os.path.exists(ppath):
			scores = dd.io.load(ppath)

		if i == 0 :
			output = scores['shap']['seq']
			shapez = output.shape
			init_beds = beds
			#raw = scores['raw']['seq']
			if 'raw' in scores:
				if one_hots is None:
					one_hots = scores['raw']['seq']
					print("one hots found")
		else:
			assert((init_beds["key"].to_numpy() == beds["key"].to_numpy()).all())
			print(scores['shap']['seq'].shape)
			assert(shapez==scores['shap']['seq'].shape)
			#assert(raw==scores['raw']['seq'])
			if 'raw' in scores:
				if one_hots is None:
					one_hots = scores['raw']['seq']
					print("one hots found")

			output += scores['shap']['seq']

		print(output.shape)
		del scores


	assert(one_hots.shape==output.shape)


	profile_scores_dict = {
			'raw': {'seq': one_hots.astype(np.int8)},
			'shap': {'seq': output.astype(np.float16)/5},
			'projected_shap': {'seq': one_hots.astype(np.int8)*(output.astype(np.float16)/5)}
			}


	os.makedirs("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/"+ddtpe+"/"+cell_type+"/merged_background/", exist_ok=True)
	output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/"+ddtpe+"/"+cell_type+"/merged_background/"+cell_type+"_folds_merged"
	dd.io.save(output_prefix+"."+itype+"_scores_new_compressed.h5",
						profile_scores_dict,
						compression='blosc')
	output_prefix_bed="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/"+ddtpe+"/"+cell_type+"/merged_background/"+cell_type+"_folds_merged"						
	beds.to_csv(output_prefix+"."+itype+"_scores_new_compressed.bed", sep="\t", header=False, index=False, compression='gzip')
							
