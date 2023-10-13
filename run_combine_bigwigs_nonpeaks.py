import pandas as pd
import os

data = pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",header=None)
ndata = data[data[1]=="GM12878"].reset_index()

#merged.GM12878_100M.profile.bw
#full_GM12878_100M.profile.bigwig
#GM12878_100M_wo_bias.bw
#full_GM12878_100M_wo_bias.bw
#{}.profile.bigwig
#/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/GM12878/GM12878_07.07.2022_bias_128_4_2356_0.5_fold_1_data_type_DNASE_PE/background_interpret/GM12878.profile.bw
#GM12878.interpreted_regions_v2.bed
all_paths = []
for i,r in ndata.iterrows():
	print(i,r[2])

	ppath = os.path.join(r[2],"background_interpret/GM12878.profile.bw")
	ppath1 = os.path.join(r[2],"background_interpret/GM12878.interpreted_regions_v2.bed")

	if os.path.isfile(ppath):
		all_paths.append(ppath)
		os.system("head "+ppath1)
	else:
		print("404 not found")

print(len(all_paths))
print(all_paths)
if len(all_paths)==5:
	output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/GM12878/merge_folds/GM12878"
	#outfile = "{}.profile.bigwig".format(output_prefix)
	outfile = "{}.profile.background.bw".format(output_prefix)

	if not os.path.isfile(outfile):
		command = "bash combine_bigwigs.sh "+all_paths[0]+" "+all_paths[1]+" "+all_paths[2]+" "+all_paths[3]+" "+all_paths[4]+" "+outfile
		print(command)
		os.system(command)


