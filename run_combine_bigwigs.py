import pandas as pd
import os

data = pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_subsample_atac.csv",header=None)
ndata = data[data[1]=="100M"].reset_index()

#merged.GM12878_100M.counts.bw
#full_GM12878_100M.counts.bigwig
#GM12878_100M_wo_bias.bw
#full_GM12878_100M_wo_bias.bw
#{}.counts.bigwig

all_paths = []
for i,r in ndata.iterrows():
	print(i,r[3])

	ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M.counts.bigwig")
	if os.path.isfile(ppath):
		all_paths.append(ppath)
	elif os.path.isfile(os.path.join(r[3],"interpret/merged.GM12878_100M.counts.bw")):
		all_paths.append(os.path.join(r[3],"interpret/merged.GM12878_100M.counts.bw"))
	elif os.path.isfile(os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M.counts.bigwig")):
		ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M.counts.bigwig")
		all_paths.append(ppath)
	else:
		print("404 not found")

print(len(all_paths))
print(all_paths)

if len(all_paths)==5:
	output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/merge_folds/full_GM12878_100M"
	#outfile = "{}.counts.bigwig".format(output_prefix)
	outfile = "{}.counts.bigwig".format(output_prefix)

	if not os.path.isfile(outfile):
		command = "bash combine_bigwigs.sh "+all_paths[0]+" "+all_paths[1]+" "+all_paths[2]+" "+all_paths[3]+" "+all_paths[4]+" "+outfile
		print(command)
		os.system(command)



all_paths = []
for i,r in ndata.iterrows():
	print(i,r[3])

	ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M.profile.bigwig")
	if os.path.isfile(ppath):
		all_paths.append(ppath)
	elif os.path.isfile(os.path.join(r[3],"interpret/merged.GM12878_100M.profile.bw")):
		all_paths.append(os.path.join(r[3],"interpret/merged.GM12878_100M.profile.bw"))
	elif os.path.isfile(os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M.profile.bigwig")):
		ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M.profile.bigwig")
		all_paths.append(ppath)
	else:
		print("404 not found")

print(all_paths)
print(len(all_paths))
if len(all_paths)==5:
	output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/merge_folds/full_GM12878_100M"
	#outfile = "{}.profile.bigwig".format(output_prefix)
	outfile = "{}.profile.bigwig".format(output_prefix)

	if not os.path.isfile(outfile):
		command = "bash combine_bigwigs.sh "+all_paths[0]+" "+all_paths[1]+" "+all_paths[2]+" "+all_paths[3]+" "+all_paths[4]+" "+outfile
		print(command)
		os.system(command)



all_paths = []
for i,r in ndata.iterrows():
	print(i,r[3])

	ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M_wo_bias.bw")
	if os.path.isfile(ppath):
		all_paths.append(ppath)
	elif os.path.isfile(os.path.join(r[3],"interpret/GM12878_100M_wo_bias.bw")):
		all_paths.append(os.path.join(r[3],"interpret/GM12878_100M_wo_bias.bw"))
	elif os.path.isfile(os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M_wo_bias.bw")):
		ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M_wo_bias.bw")
		all_paths.append(ppath)
	else:
		print("404 not found")

print(all_paths)
print(len(all_paths))
if len(all_paths)==5:
	output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/merge_folds/full_GM12878_100M"
	#outfile = "{}.counts.bigwig".format(output_prefix)
	outfile = "{}_wo_bias.bigwig".format(output_prefix)

	if not os.path.isfile(outfile):
		command = "bash combine_bigwigs.sh "+all_paths[0]+" "+all_paths[1]+" "+all_paths[2]+" "+all_paths[3]+" "+all_paths[4]+" "+outfile
		print(command)
		os.system(command)


all_paths = []
for i,r in ndata.iterrows():
	print(i,r[3])

	ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M_w_bias.bw")
	if os.path.isfile(ppath):
		all_paths.append(ppath)
	elif os.path.isfile(os.path.join(r[3],"interpret/GM12878_100M_w_bias.bw")):
		all_paths.append(os.path.join(r[3],"interpret/GM12878_100M_w_bias.bw"))
	elif os.path.isfile(os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M_w_bias.bw")):
		ppath = os.path.join(r[3],"chrombpnet_model/interpret/full_GM12878_100M_w_bias.bw")
		all_paths.append(ppath)
	else:
		print("404 not found")

print(all_paths)
print(len(all_paths))
if len(all_paths)==5:

	output_prefix="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/GM12878_100M/merge_folds/full_GM12878_100M"
	outfile = "{}_w_bias.bigwig".format(output_prefix)

	if not os.path.isfile(outfile):
		command = "bash combine_bigwigs.sh "+all_paths[0]+" "+all_paths[1]+" "+all_paths[2]+" "+all_paths[3]+" "+all_paths[4]+" "+outfile
		print(command)
		os.system(command)



