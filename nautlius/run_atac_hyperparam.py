import os
import shutil
from os.path import exists


#bias_filt,dil_layers,seed,bias_th
hyper_params=[["128","4","1234", "0.5"], ["128","4","2356","0.5"], ["128","4","5678","0.5"]]
#hyper_params.extend([["256","5","1234","0.5"], ["256","5","2356","0.5"], ["256","5","5678","0.5"]])
#hyper_params.extend([["512","8","1234","0.5"], ["512","8","2356","0.5"], ["512","8","5678","0.5"]])


job_meta = open("job.yml").read()

for param in hyper_params:

	with open("job_scripts/run_"+param[0]+"_"+param[1]+"_"+param[2]+"_"+param[3]+".yml", "w") as f:
		f.write(job_meta.replace("$experiment","K562").replace("$bias_filt",param[0]).replace("$bias_dil",param[1]).replace("$seed",param[2]).replace("$bias_th",param[3]))
		command = "kubectl create -f "+"job_scripts/run_"+param[0]+"_"+param[1]+"_"+param[2]+"_"+param[3]+".yml"
	
	print(command)
	os.system(command)

