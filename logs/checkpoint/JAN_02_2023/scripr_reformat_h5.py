
import pandas as pd
import os

data = pd.read_csv("model_dir_atac.csv",header=None)
ddtpe="ATAC"
ddtpen=ddtpe+"_PE"
cell_types=["HEPG2", "K562", "GM12878", "H1ESC", "IMR90"]
#cell_types=[ "H1ESC", "IMR90"]
#cell_types=["H1ESC", "IMR90"]

for cell_type in cell_types:
	ndata = data[data[1]==cell_type].reset_index()
	output_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+cell_type+"/interpret_upload/"
	if not os.path.isdir(output_dir):
		os.mkdir(output_dir)

	for i,r in ndata.iterrows():
		#print(i,r[2])

		output_dir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+cell_type+"/interpret_upload/"+r[0]
		if not os.path.isdir(output_dir):
			os.mkdir(output_dir)


		itype="counts"
		output_prefix=os.path.join(r[2],"chrombpnet_model/interpret_all/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		output_prefix_bed=os.path.join(r[2],"chrombpnet_model/interpret_all/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		outputfile=output_dir+"/"+cell_type+"_"+itype

		if not os.path.isfile(outputfile+"_attribs_reformatted.h5"):
			command = "python reformat_h5s.py -h5py1 "+output_prefix+" -r1 "+output_prefix_bed+" -o "+outputfile
			print(command)
			os.system(command)

		itype="profile"
		output_prefix=os.path.join(r[2],"chrombpnet_model/interpret_all/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		output_prefix_bed=os.path.join(r[2],"chrombpnet_model/interpret_all/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		outputfile=output_dir+"/"+cell_type+"_"+itype


		if not os.path.isfile(outputfile+"_attribs_reformatted.h5"):
			command = "python reformat_h5s.py -h5py1 "+output_prefix+" -r1 "+output_prefix_bed+" -o "+outputfile
			print(command)
			os.system(command)

	#	break
	#break
