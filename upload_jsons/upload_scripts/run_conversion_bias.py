import pandas as pd
import os

#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",",header=None)
model_atac = pd.read_csv("bias_models_atlas.csv",sep=",",header=None)
#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",",header=None)


for i,r in model_atac.iterrows():
	fold = r[0]	
	name = r[1]
	model_path = r[2]
	input_path=os.path.join(model_path,"bias_model/bias.h5")
	output_path=os.path.join(model_path,"bias_model/new_model_formats_v1/bias")
	output_dir=os.path.join(model_path,"bias_model/new_model_formats_v1/")
	file_path="bias"
	if not os.path.isfile(input_path):
		print("ERROR bias model not found")
		print(input_path)
		continue

	print(output_path)
	if not os.path.isfile(output_path+".tar"):

		os.makedirs(os.path.join(model_path,"bias_model/new_model_formats_v1/"), exist_ok=True)
		command = "CUDA_VISIBLE_DEVICES=2 python get_new_tf_model_format.py -i "+input_path+" -o "+output_dir+" -f "+file_path
		print(command)	
		os.system(command)

