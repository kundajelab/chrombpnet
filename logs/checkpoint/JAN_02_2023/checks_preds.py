import h5py
import pandas as pd

ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/preds_upload/average_preds_with_ccre_vf/HEPG2.mean_preds_w_bias_predictions.h5"
print(ifile)
data = h5py.File(ifile, "r")
print(data["predictions"]["logits"].shape)
print(data["predictions"]["logcounts"].shape)


print(data["coords"]["coords_chrom"][0], data["coords"]["coords_start_dset"][0], data["coords"]["coords_end_dset"][0])
print(data["coords"]["coords_chrom"][-1], data["coords"]["coords_start_dset"][-1], data["coords"]["coords_end_dset"][-1])

ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/preds_upload/average_preds_with_ccre_vf/HEPG2.mean_preds_wo_bias_predictions.h5"
print(ifile)
data = h5py.File(ifile, "r")
print(data["predictions"]["logits"].shape)
print(data["predictions"]["logcounts"].shape)


print(data["coords"]["coords_chrom"][0], data["coords"]["coords_start_dset"][0], data["coords"]["coords_end_dset"][0])
print(data["coords"]["coords_chrom"][-1], data["coords"]["coords_start_dset"][-1], data["coords"]["coords_end_dset"][-1])


#data = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/preds_upload/fold_0/HEPG2_w_bias_all_regions.bed", sep='\t', header=None)
data = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/preds_upload/average_preds_with_ccre_vf/filtered.regions.bed.gz", sep='\t', header=None)
print(data.shape)
print(data.head(2))
print(data.tail(2))

data = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/preds_upload/average_preds_with_ccre_vf/input.regions.bed.gz", sep='\t', header=None)
#data = data[~(data[0]=="chrM")]
print(data.shape)
print(data.head(2))
print(data.tail(2))

print(data.iloc[0,0], data.iloc[0,1]+data.iloc[0,9]-500,  data.iloc[0,1]+data.iloc[0,9]-500+1000)
print(data.iloc[-1,0], data.iloc[-1,1]+data.iloc[-1,9]-500,  data.iloc[-1,1]+data.iloc[-1,9]-500+1000)


data = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/interpret_upload/average_preds/mean_folds.inputs.bed.gz", sep="\t", header=None)
print(data.shape)
print(data.head(2))
print(data.tail(2))




for i in range(0,5):
	ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/preds_upload/fold_{}/HEPG2_w_bias_all_with_ccre_predictions.h5".format(str(i))
	print(ifile)
	data = h5py.File(ifile, "r")
	print(data["predictions"]["logits"].shape)
	print(data["predictions"]["logcounts"].shape)

for i in range(0,5):
	ifile="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/HEPG2/preds_upload/fold_{}/HEPG2_wo_bias_all_with_ccre_predictions.h5".format(str(i))
	print(ifile)
	data = h5py.File(ifile, "r")
	print(data["predictions"]["logits"].shape)
	print(data["predictions"]["logcounts"].shape)



