import pandas as pd
import numpy as np
from tqdm import tqdm
import deepdish 
import argparse



NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


parser = argparse.ArgumentParser(description="subsample scores")
parser.add_argument("-c", "--celltype", type=str, required=True)
#parser.add_argument("-mn", "--modelname", type=str, required=True)
#parser.add_argument("-dt", "--datatype", type=str, required=True)
parser.add_argument("-md", "--modeldir", type=str, default="", required=False)
args = parser.parse_args()

celltype=args.celltype
modelname="chrombpnet_model"
#data_type=args.datatype
modeldir=args.modeldir

#print(celltype,modelname,data_type,modeldir)
print(celltype,modeldir)
#celltype="K562"
#modelname="GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0"
#data_type="ATAC_PE"
#modeldir="nautilus_runs"

#/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/GM12878_100M_07.19.2022_bias_transfer_1234_fold_0_data_type_ATAC_PE/chrombpnet_model/interpret/full_GM12878_100M.interpreted_regions_counts.bed



merged = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC"+"/"+celltype+"/"+modeldir+"/"+modelname+"/interpret/full_"+celltype+".interpreted_regions_counts.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
bed_of_interest = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/GM12878/data/peaks_no_blacklist.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
counts_full_h5py_file="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC"+"/"+celltype+"/"+modeldir+"/"+modelname+"/interpret/full_"+celltype+".counts_scores.h5"
profile_full_h5py_file="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC"+"/"+celltype+"/"+modeldir+"/"+modelname+"/interpret/full_"+celltype+".profile_scores.h5"
#output_prefix="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/"+data_type+"/"+celltype+"/"+modelname+"/SIGNAL/"+celltype+"_full"
output_prefix="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC"+"/"+celltype+"/"+celltype+"_full"
print(output_prefix)

boi = bed_of_interest[["chr", "start", "end", "summit"]].to_numpy().tolist()
merged_val = merged[["chr", "start", "end", "summit"]].to_numpy().tolist()

indices=[]
for i, val in enumerate(tqdm(merged_val)):
	if val in boi:
		indices.append(i)

print(len(indices))
print(len(merged_val))
print(len(boi))

merged.iloc[indices].to_csv("{}.interpreted_regions.bed".format(output_prefix), header=False, index=False, sep="\t")
# subsample counts

scores = deepdish.io.load(counts_full_h5py_file)
sub_scores = {
            'shap': {'seq': scores['shap']['seq'][indices]},
        }

print(sub_scores['shap']['seq'].shape)
deepdish.io.save("{}.counts_scores.h5".format(output_prefix),
                    sub_scores,
                    compression='blosc')


del sub_scores, scores

# subsample profile

scores = deepdish.io.load(profile_full_h5py_file)
sub_scores = {
            'shap': {'seq': scores['shap']['seq'][indices]},
        }

print(sub_scores['shap']['seq'].shape)
deepdish.io.save("{}.profile_scores.h5".format(output_prefix),
                    sub_scores,
                    compression='blosc')


del sub_scores, scores