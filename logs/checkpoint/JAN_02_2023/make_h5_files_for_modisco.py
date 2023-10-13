import pandas as pd
import numpy as np
from tqdm import tqdm
import deepdish 
import argparse



NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]


#parser = argparse.ArgumentParser(description="subsample scores")
#parser.add_argument("-c", "--celltype", type=str, required=True)
#parser.add_argument("-mn", "--modelname", type=str, required=True)
#parser.add_argument("-dt", "--datatype", type=str, required=True)
#parser.add_argument("-md", "--modeldir", type=str, default="", required=False)
#args = parser.parse_args()

#celltype=args.celltype
#modelname=args.modelname
#data_type=args.datatype
#modeldir=args.modeldir

#print(celltype,modelname,data_type,modeldir)
#celltype="GM12878"
#modelname="GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0"
#data_type="ATAC_PE"
#modeldir="nautilus_runs"

#merged = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"+data_type+"/"+celltype+"/"+modeldir+"/"+modelname+"/interpret/merged."+celltype+".interpreted_regions.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
#bed_of_interest = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"+data_type+"/"+celltype+"/data/peaks_no_blacklist.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
#counts_full_h5py_file="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"+data_type+"/"+celltype+"/"+modeldir+"/"+modelname+"/interpret/merged."+celltype+".counts_scores.h5"
#profile_full_h5py_file="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"+data_type+"/"+celltype+"/"+modeldir+"/"+modelname+"/interpret/merged."+celltype+".profile_scores.h5"
#output_prefix="/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/"+data_type+"/"+celltype+"/"+modelname+"/SIGNAL/"+celltype+"_full"
#output_prefix="/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/modisco/"+data_type.split("_")[0]+"/"+celltype+"/"+celltype+"_full"

data_type="ATAC_PE"

celltype="GM12878"
#celltype="IMR90"
#merged = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_07.07.2022_bias_128_4_2356_0.5_fold_1_data_type_ATAC_PE/chrombpnet_model/interpret/full_GM12878.interpreted_regions_counts.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
#merged = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_07.07.2022_bias_128_4_1234_0.5_fold_2_data_type_ATAC_PE/chrombpnet_model/interpret/full_GM12878.interpreted_regions_counts.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)


merged = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878/GM12878_07.08.2022_bias_128_4_1234_0.4_fold_1_data_type_ATAC_PE/chrombpnet_model/interpret/full_GM12878.interpreted_regions_counts.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
#merged = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/IMR90/IMR90_07.17.2022_bias_128_4_1234_0.3_fold_1_data_type_ATAC_PE/chrombpnet_model/interpret/full_IMR90.interpreted_regions_counts.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
bed_of_interest = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"+data_type+"/"+celltype+"/data/peaks_no_blacklist.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)


counts_full_h5py_file="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+celltype+"/merge_folds/"+celltype+"_folds_merged.counts_scores.h5"
profile_full_h5py_file="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+celltype+"/merge_folds/"+celltype+"_folds_merged.profile_scores.h5"
output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+celltype+"/merge_folds/in_peaks"


print(output_prefix)


boi = bed_of_interest[["chr", "start", "end", "summit"]].to_numpy().tolist()
merged_val = merged[["chr", "start", "end", "summit"]].to_numpy().tolist()

indices=[]
for i, val in enumerate(tqdm(merged_val)):
	if val in boi:
		indices.append(i)

print(len(indices))
print(len(merged_val))

merged.iloc[indices].to_csv("{}.interpreted_regions.bed".format(output_prefix), header=False, index=False, sep="\t")
# subsample counts

scores = deepdish.io.load(counts_full_h5py_file)
sub_scores = {
            'raw': {'seq': scores['raw']['seq'][indices]},
            'shap': {'seq': scores['shap']['seq'][indices]},
            'projected_shap': {'seq': scores['projected_shap']['seq'][indices]}
        }

print(sub_scores['raw']['seq'].shape)
deepdish.io.save("{}.counts_scores.h5".format(output_prefix),
                    sub_scores,
                    compression='blosc')


del sub_scores, scores

# subsample profile

scores = deepdish.io.load(profile_full_h5py_file)
sub_scores = {
            'raw': {'seq': scores['raw']['seq'][indices]},
            'shap': {'seq': scores['shap']['seq'][indices]},
            'projected_shap': {'seq': scores['projected_shap']['seq'][indices]}
        }

print(sub_scores['raw']['seq'].shape)
deepdish.io.save("{}.profile_scores.h5".format(output_prefix),
                    sub_scores,
                    compression='blosc')


del sub_scores, scores
