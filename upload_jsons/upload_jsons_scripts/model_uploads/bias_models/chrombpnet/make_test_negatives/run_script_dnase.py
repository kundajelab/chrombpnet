import pandas as pd
import os

#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_scripts/model_dir_dnase_v2.1_bias.csv",sep=",", header=None)
model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/v1/model_dir_dnase_v2.1.csv",sep=",", header=None)


print(model_atac.head())

for i,r in model_atac.iterrows():

	if r[1] != "H1ESC":
		continue
	if r[1] in ["HEPG2", "K562"]:
		tag="DNASE_PE"
		mdir=r[1]
	else:
		tag="DNASE_SE"
	#print(r)
	if os.path.isfile(os.path.join(r[2], "train_test_regions_may_7_2024/nonpeaks.testset.bed.gz")):
		try:
			tdata = pd.read_csv(os.path.join(r[2], "train_test_regions_may_7_2024/nonpeaks.testset.bed.gz"))
			continue
		except:
			pass

	print(os.path.join(r[2], "train_test_regions_may_7_2024/nonpeaks.testset.bed.gz"))
	if r[0].split("_")[-1] == "0":
		command = "bash script_dnase.sh "+r[2]+" "+r[1]+" "+r[0]+" "+tag
	else:
		command = "bash script_dnase.sh "+r[2]+" "+r[1]+" "+r[0]+" "+tag+" "+"_"+str(r[0].split("_")[-1])

	print(command)
	os.system(command)	

#/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/H1ESC/n
