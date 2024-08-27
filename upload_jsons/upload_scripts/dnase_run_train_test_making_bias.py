import pandas as pd
import os

#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",", header=None)
#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",", header=None)
#model_atac=pd.read_csv("bias_models_atlas.csv", sep=',', header=None)
model_atac=pd.read_csv("model_dir_dnase_v2.1_bias.csv", sep=',', header=None)




encode_id = {"K562": "ENCSR868FGK",
        "GM12878": "ENCSR637XSC",
        "HEPG2": "ENCSR291GJU",
        "IMR90": "ENCSR200OML",
        "H1ESC": "ENCDUMMY"}

encode_id = {"K562": "ENCSR000EOT",
        "GM12878": "ENCSR000EMT",
        "HEPG2": "ENCSR149XIL",
        "IMR90": "ENCSR477RTP",
        "H1ESC": "ENCSR000EMU"}


for i,r in model_atac.iterrows():
	fold = r[0]
	name = r[1]
	model_path = r[2]

	#input_peaks=os.path.join(model_path,"chrombpnet_model/filtered.peaks.bed")
	input_nonpeaks=os.path.join(model_path,"bias_model/filtered.bias_nonpeaks.bed")
	#test_nonpeaks="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"+encode_id[name]+"/negatives_data/test/test."+fold+".filtered.negatives_with_summit.bed"
	fold="/mnt/lab_data2/anusri/chrombpnet/splits/"+fold+".json"
	output_path=os.path.join(model_path,"train_test_regions_bias_may_7_2024/")

	if not os.path.isfile(input_nonpeaks):
		cellline=input_nonpeaks.split("/")[10]
		biasth=input_nonpeaks.split("/")[11].split("_")[6]
		foldn=input_nonpeaks.split("/")[11].split("_")[8]
		#print(cellline,biasth,foldn)
		if cellline in ["K562", "HEPG2"]:
			ddatype="DNASE_PE"
		elif cellline in ["H1ESC"]:
			ddatype="DNASE_SE"
		else:
			print(cellline)
			break
		outputdir=os.path.join(model_path,"bias_model/newgen/")
		if not  os.path.isfile(os.path.join(model_path,"bias_model/newgen/filtered.bias_nonpeaks.bed")):
			os.makedirs(outputdir, exist_ok=True)
			print(outputdir)
			command = "bash make_missing_bed_regions.sh "+cellline+" "+biasth+" "+foldn+" "+outputdir+" "+ddatype
			os.system(command)
			print(command)
			input_nonpeaks=os.path.join(model_path,"bias_model/newgen/filtered.bias_nonpeaks.bed")
		else:
			input_nonpeaks=os.path.join(model_path,"bias_model/newgen/filtered.bias_nonpeaks.bed")

	if not os.path.isfile(output_path+"nonpeaks.validationset.bed.gz"):
		print(output_path)
		os.makedirs(output_path, exist_ok=True)
		command=["python get_train_test_regions_bias.py "]+["-inp"]+[input_nonpeaks]+["-f"]+[fold]+["-o"]+[output_path]
		command = " ".join(command)
		print(command)
		os.system(command)
	else:
		print(output_path)
