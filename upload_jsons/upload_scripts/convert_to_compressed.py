import pandas as pd
import os

for name in [ "HEPG2", "IMR90", "H1ESC"]:

	input_path="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/data/peaks_no_blacklist.bed"
	if os.path.isfile(input_path):
		data = pd.read_csv(input_path, sep="\t", header=None)
		data.to_csv(input_path+".gz", compression="gzip", sep="\t", index=False, header=False)

	input_path="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/negatives_data/negatives_with_summit.bed"
	if os.path.isfile(input_path):
		data = pd.read_csv(input_path, sep="\t", header=None)
		data.to_csv(input_path+".gz", compression="gzip", sep="\t", index=False, header=False)


	input_path="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/negatives_data_1/negatives_with_summit.bed"
	if os.path.isfile(input_path):
		data = pd.read_csv(input_path, sep="\t", header=None)
		data.to_csv(input_path+".gz", compression="gzip", sep="\t", index=False, header=False)


	input_path="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/negatives_data_2/negatives_with_summit.bed"
	if os.path.isfile(input_path):
		data = pd.read_csv(input_path, sep="\t", header=None)
		data.to_csv(input_path+".gz", compression="gzip", sep="\t", index=False, header=False)


	input_path="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/negatives_data_3/negatives_with_summit.bed"
	if os.path.isfile(input_path):
		data = pd.read_csv(input_path, sep="\t", header=None)
		data.to_csv(input_path+".gz", compression="gzip", sep="\t", index=False, header=False)


	input_path="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/negatives_data_4/negatives_with_summit.bed"
	if os.path.isfile(input_path):
		data = pd.read_csv(input_path, sep="\t", header=None)
		data.to_csv(input_path+".gz", compression="gzip", sep="\t", index=False, header=False)

