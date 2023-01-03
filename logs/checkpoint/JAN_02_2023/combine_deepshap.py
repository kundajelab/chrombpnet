import pandas as pd
import os
import deepdish as dd

data = pd.read_csv("model_dir_atac.csv",header=None)
ndata = data[data[1]=="K562"].reset_index()

for i,r in ndata.iterrows():
	print(i,r[2])
	ppath = os.path.join(r[2],"chrombpnet_model/interpret/full_K562.counts_scores.h5")
	if os.path.exists(ppath):
		scores = dd.io.load(ppath)
	else:
		ppath = os.path.join(r[2],"interpret/merged.K562.counts_scores.h5")
		scores = dd.io.load(ppath)


	if i == 0 :
		output = scores['shap']['seq']
		shapez = output.shape
	else:
		print(scores['shap']['seq'].shape)
		assert(shapez==scores['shap']['seq'].shape)
		output += scores['shap']['seq']

	print(output.shape)
	del scores


for i,r in ndata.iterrows():
	print(i,r[2])
	ppath = os.path.join(r[2],"chrombpnet_model/interpret/full_K562.counts_scores.h5")
	if os.path.exists(ppath):
		scores = dd.io.load(ppath)
	else:
		ppath = os.path.join(r[2],"interpret/merged.K562.counts_scores.h5")
		scores = dd.io.load(ppath)
	break

counts_scores_dict = {
            'raw': {'seq': scores['shap']['raw']},
            'shap': {'seq': output},
            'projected_shap': {'seq': scores['shap']['raw']*output}
        }


output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/K562/merge_folds/K562_folds_merged"
dd.io.save("{}.counts_scores.h5".format(output_prefix),
                    counts_scores_dict,
                    compression='blosc')
