import h5py
import pandas as pd

ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/K562_{}_attribs_reformatted.h5".format("counts")
print(ifile)
data = h5py.File(ifile, "r")
print(data["attributions"]["shap"].shape)

print(data["coords"]["coords_chrom"][0], data["coords"]["coords_start_dset"][0], data["coords"]["coords_end_dset"][0])
print(data["coords"]["coords_chrom"][-1], data["coords"]["coords_start_dset"][-1], data["coords"]["coords_end_dset"][-1])

ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/K562_{}_attribs_reformatted.h5".format("profile")
print(ifile)
data = h5py.File(ifile, "r")
print(data["attributions"]["shap"].shape)

print(data["coords"]["coords_chrom"][0], data["coords"]["coords_start_dset"][0], data["coords"]["coords_end_dset"][0])
print(data["coords"]["coords_chrom"][-1], data["coords"]["coords_start_dset"][-1], data["coords"]["coords_end_dset"][-1])



data = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/mean_folds.inputs.bed.gz", sep="\t", header=None)
print(data.shape)
print(data.head(2))
print(data.tail(2))

print(data.iloc[0,0], data.iloc[0,1]+data.iloc[0,9]-1057,  data.iloc[0,1]+data.iloc[0,9]-1057+2114)
print(data.iloc[-1,0], data.iloc[-1,1]+data.iloc[-1,9]-1057,  data.iloc[-1,1]+data.iloc[-1,9]-1057+2114)



ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/merge_folds_new_may_05_24/in_peaks.counts_scores_new_compressed.h5"
print(ifile)
data = h5py.File(ifile, "r")
print(data["shap"]['seq'].shape)

data=pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/modisco.inputs.bed.gz", sep='\t', header=None)
print(data.shape)
print(data.head(2))
print(data.tail(2))

#print(data.[0,0], data[0][1]+data[0][9]-1057,data[0][1]+data[0][9]-1057+2114)
#print(data[-1][0], data[-1][1]+data[0][9]-1057,data[-1][1]+data[-1][9]-1057+2114)



for i in range(0,5):
	ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/fold_{}/K562_{}_attribs_reformatted.h5".format(str(i),"counts")
	print(ifile)
	data = h5py.File(ifile, "r")
	print(data["attributions"]["shap"].shape)



for i in range(0,5):
	ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/fold_{}/K562_{}_attribs_reformatted.h5".format(str(i),"profile")
	print(ifile)
	data = h5py.File(ifile, "r")
	print(data["attributions"]["shap"].shape)



data=pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/K562/interpret_upload/average_preds/per_folds.inputs.bed.gz", sep='\t', header=None)
print(data.shape)
print(data.head(2))
print(data.tail(2))

