import pandas as pd
import os

model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",",header=None)
#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",",header=None)

for i,r in model_atac.iterrows():
	fold = r[0]	
	name = r[1]
	model_path = r[2]
	input_path=os.path.join(model_path,"chrombpnet_model/chrombpnet.h5")
	output_path=os.path.join(model_path,"chrombpnet_model/new_model_formats/chrombpnet")
	print(output_path)
	if not os.path.isfile(output_path+".tar"):

		os.makedirs(os.path.join(model_path,"chrombpnet_model/new_model_formats/"), exist_ok=True)
		command = "CUDA_VISIBLE_DEVICES=5 python get_new_tf_model_format.py -i "+input_path+" -o "+output_path
		print(command)	
		os.system(command)

	input_path=os.path.join(model_path,"chrombpnet_model/chrombpnet_wo_bias.h5")
	output_path=os.path.join(model_path,"chrombpnet_model/new_model_formats/chrombpnet_wo_bias")

	if not os.path.isfile(output_path+".tar"):
		os.makedirs(os.path.join(model_path,"chrombpnet_model/new_model_formats/"), exist_ok=True)
		command = "CUDA_VISIBLE_DEVICES=5 python get_new_tf_model_format.py -i "+input_path+" -o "+output_path
		print(command)	
		os.system(command)

	input_path=os.path.join(model_path,"chrombpnet_model/bias_model_scaled.h5")
	output_path=os.path.join(model_path,"chrombpnet_model/new_model_formats/bias_model_scaled")

	if not os.path.isfile(output_path+".tar"):
		os.makedirs(os.path.join(model_path,"chrombpnet_model/new_model_formats/"), exist_ok=True)
		command = "CUDA_VISIBLE_DEVICES=5 python get_new_tf_model_format.py -i "+input_path+" -o "+output_path
		print(command)	
		os.system(command)


