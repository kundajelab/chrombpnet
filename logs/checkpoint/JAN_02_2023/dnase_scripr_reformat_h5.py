
import pandas as pd
import os

#data = pd.read_csv("model_dir_atac.csv",header=None)
#data = pd.read_csv("model_dir_dnase.csv",header=None)
data = pd.read_csv("v1/model_dir_dnase_v2_interpret.csv",header=None)

ddtpe="DNASE"
ddtpen=ddtpe+"_PE"
#cell_types=["HEPG2", "K562", "GM12878", "H1ESC", "IMR90"]
#cell_types=["HEPG2", "K562"]
cell_types=[ "H1ESC"]
#cell_types=["IMR90", "GM12878"]

  
for cell_type in cell_types:

	ndata = data[data[1]==cell_type].reset_index()
	cell_type = cell_type+"_new"
	output_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/interpret_upload/"
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)

	for i,r in ndata.iterrows():
		#print(i,r[2])

		output_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/interpret_upload/"
		if not os.path.isdir(output_dir):
			os.mkdir(output_dir)

		output_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/interpret_upload/"+r[0]
		if not os.path.isdir(output_dir):
			os.mkdir(output_dir)

		itype="counts"

		#rpath="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/"+r[0]
		rpath="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/"+r[2].split("/")[-2]
		output_prefix=os.path.join(rpath,"chrombpnet_model/interpret_all_with_ccre/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		output_prefix_bed=os.path.join(rpath,"chrombpnet_model/interpret_all_with_ccre/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		#output_prefix=os.path.join(rpath,"interpret_all_with_ccre/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		#output_prefix_bed=os.path.join(rpath,"interpret_all_with_ccre/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		outputfile=output_dir+"/"+cell_type+"_"+itype

		command = "python reformat_h5s.py -h5py1 "+output_prefix+" -r1 "+output_prefix_bed+" -o "+outputfile
		print(command)
		os.system(command)

		itype="profile"

		#rpath="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/"+r[2].split("/")[-2]
		output_prefix=os.path.join(rpath,"chrombpnet_model/interpret_all_with_ccre/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		output_prefix_bed=os.path.join(rpath,"chrombpnet_model/interpret_all_with_ccre/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		#rpath="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/"+r[0]
		#output_prefix=os.path.join(rpath,"interpret_all_with_ccre/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		#output_prefix_bed=os.path.join(rpath,"interpret_all_with_ccre/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		outputfile=output_dir+"/"+cell_type+"_"+itype

		command = "python reformat_h5s.py -h5py1 "+output_prefix+" -r1 "+output_prefix_bed+" -o "+outputfile
		print(command)
		os.system(command)

	#	break
	#break
