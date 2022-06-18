import os
import shutil
from os.path import exists


#bias_filt,dil_layers,seed,bias_th

#params=[["128","4","1234", "0.5", "K562", "fold_0"], ["128","4","1234", "0.5", "K562", "fold_1"], ["128","4","1234", "0.5", "K562", "fold_2"], ["128","4","1234", "0.5", "K562", "fold_3"], ["128","4","1234", "0.5", "K562", "fold_4"]]
params=[["128","4","1234", "0.5", "K562", "fold_0"]]
hyper_params=[]
for cell_line in ["HEPG2", "IMR90", "H1ESC", "GM12878", "K562"]:
	for param in params:
		if cell_line!="IMR90":
			#continue
			hyper_params.append(param[0:4]+[cell_line]+[param[5]])
		else:
			hyper_params.append(param[0:3]+["0.4"]+[cell_line]+[param[5]])

print(hyper_params)
job_meta = open("job_dataset_atac.yml").read()

for param in hyper_params:

	with open("job_scripts_apr_8/run_atac_stride_50_"+param[4]+"_"+param[5]+"_"+param[0]+"_"+param[1]+"_"+param[2]+"_"+param[3]+".yml", "w") as f:
		f.write(job_meta.replace("$cell_line",param[4]).replace("$dataset",param[4].lower()).replace("$fold_num_name",param[5].replace("_",".")).replace("$fold_num",param[5]).replace("$bias_filt",param[0]).replace("$bias_dil",param[1]).replace("$seed",param[2]).replace("$bias_th",param[3]))
		command = "kubectl create -f "+"job_scripts_apr_8/run_atac_stride_50_"+param[4]+"_"+param[5]+"_"+param[0]+"_"+param[1]+"_"+param[2]+"_"+param[3]+".yml"
	
	print(command)
	os.system(command)

