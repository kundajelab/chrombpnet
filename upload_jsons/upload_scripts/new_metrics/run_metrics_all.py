
import pandas as pd
import os

model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",", header=None)
#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",")

encode_id = {"K562": "ENCSR868FGK",
	"GM12878": "ENCSR637XSC",
	"HEPG2": "ENCSR291GJU",
	"IMR90": "ENCSR200OML"}


for i,r in model_atac.iterrows():
	fold = r[0]
	name = r[1]
	model_path = r[2]
	if name != "IMR90":
		continue

	bigwig_path="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/data/"+name+"_unstranded.bw"
	nonpeaks="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/"+encode_id[name]+"/negatives_data/test/test."+fold+".filtered.negatives_with_summit.bed"
	fold="/mnt/lab_data2/anusri/chrombpnet/splits/"+fold+".json"

	chrombpnet=os.path.join(model_path,"chrombpnet_model/chrombpnet.h5")
	chrombpnet_nb=os.path.join(model_path,"chrombpnet_model/chrombpnet_wo_bias.h5")
	if os.path.isfile(chrombpnet_nb):
		if os.path.isfile(os.path.join(model_path, "chrombpnet_model/new_test_metrics/chrombpnet_metrics.json")):
			continue
		print(os.path.join(model_path, "new_test_metrics/chrombpnet_metrics.json"))
		command="CUDA_VISIBLE_DEVICES=1 bash new_metrics.sh "+bigwig_path+" "+nonpeaks+" "+fold+" "+os.path.join(model_path,"chrombpnet_model/")
		print(command)
		os.system(command)
