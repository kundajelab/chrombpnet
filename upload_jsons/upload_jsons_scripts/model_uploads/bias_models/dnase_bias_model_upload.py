import os
import dnase_bias_upload_utils as upload_utils
import json
import pandas as pd

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"
#output_dir =  "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/DNASE/stage1/jul_17_2023/"
output_dir = "dnase_production_uploads/"

encids = os.listdir(odir)
#encids = open("../chromatin_atlas_atac/test_encid.txt").readlines()
#encids = [line.strip() for line in encids]

model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_scripts/bias_models_atlas.csv",sep=",", header=None)

encode_id = {"HEPG2": "ENCSR149XIL",
        "ENCSR283TME": "ENCSR283TME",
        "ENCSR880CUB": "ENCSR880CUB",
        "ENCSR146KFX": "ENCSR146KFX"}
        
data_to_bam = {"HEPG2": ["ENCFF474LSZ", "ENCFF839SPF"],
        "ENCSR283TME": ["ENCFF906TID"],
        "ENCSR880CUB": ["ENCFF966VCS"],
        "ENCSR146KFX": ["ENCFF639ETA"]}
        

def main_fetch_preprocessing_files(encid, args_json, bam_ids, name):
	
	success_flag = False
	
	args_json["upload bias"] = True
	args_json["bias model encid"] = encid 

	# find the bams input
	preprocessing_path = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/"+name+"/data/"+name+"_unstranded.bw"
	preprocessing_path_oak = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"+encid+"/preprocessing/bigWigs/"+encid+".bigWig"
	if os.path.isfile(preprocessing_path):	
		args_json["experiment"] = encid
		args_json["bam files"] = bam_ids
		args_json["assay"] = "DNase-seq"
		args_json["observed signal profile bigWig"] = preprocessing_path
		success = True
	elif os.path.isfile(preprocessing_path_oak):
		args_json["experiment"] = encid
		args_json["bam files"] = bam_ids
		args_json["assay"] = "DNase-seq"
		args_json["observed signal profile bigWig"] = preprocessing_path_oak
		success = True		
	else:
		success = False
	
	return	success, args_json

def main_fetch_bias_model_files(encid, args_json, models_path):
	success = False
	args_json["bias models tar"] = {}
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/bias.models.README"
	assert(os.path.isfile(readme_file))
	args_json["bias models tar"]["file.paths"] = [(readme_file, "README.md")]
	#args_json["bias models tar"]["logs.bias.models."+encid] = {"file.paths": None}

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_bias_models(odir, models_path[i], encid, i)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["bias models tar"]["fold_"+str(i)] = {}
		args_json["bias models tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["bias models tar"]["fold_"+str(i)]["logs.bias.models.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		#print(log_paths)
		# 9 log file expected per model
		assert(len(log_paths) == 7)	
		assert(len(data_paths) == 2)		
	success=True
	return success, args_json
	

def main_fetch_bias_training_files(encid, args_json, models_path, name):
	success = False
	
	# find the training test regions
	args_json["bias training and test regions tar"] = {}	
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/bias.training_test_regions.README"
	assert(os.path.isfile(readme_file))
	args_json["bias training and test regions tar"]["file.paths"] = [(readme_file, "README.md")]

	if name == "HEPG2":
		main_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_PE/"
		input_peaks = os.path.join(main_dir, name + "/data/peaks_no_blacklist.bed.gz")
	else:
		main_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/DNASE_SE/"
		input_peaks = os.path.join(odir, encid + "/preprocessing/downloads/peaks.bed.gz")

	print(input_peaks)
	if os.path.isfile(input_peaks):
		args_json["bias training and test regions tar"]["file.paths"].append((input_peaks,"peaks.all_input_regions."+encid+".bed.gz"))		
	else:
		success = False
		return success, args_json
		
	# log files preprocessing and peak-calling
	if name == "HEPG2":
		log_paths = upload_utils.bias_fetch_preprocessing_log_files_set_1(odir, encid, main_dir, name)
		assert(len(log_paths) == 6)		
	else:
		log_paths = upload_utils.bias_fetch_preprocessing_log_files_set_2(odir, encid, main_dir, name)
		assert(len(log_paths) == 8)	

		
	args_json["bias training and test regions tar"]["logs.bias.training_test_regions."+encid] = {"file.paths": log_paths}
	

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_training_data_bias(odir, models_path[i], encid, i, main_dir, name)
		#print(data_paths)
		args_json["bias training and test regions tar"]["fold_"+str(i)] = {}
		args_json["bias training and test regions tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["bias training and test regions tar"]["fold_"+str(i)]["logs.bias.training_test_regions.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		#print(log_paths)
		#print(log_paths)
		print(data_paths)
		assert(len(data_paths) == 5)
		assert(len(log_paths) == 2)	

		#if len(data_paths) != 3:
		#	success = False
		#	return success, args_json
	
	success = True
	return success, args_json

	
		
if __name__ == "__main__":

	# define readmes specfic to bias model
	for name in ["HEPG2", "ENCSR283TME", "ENCSR880CUB", "ENCSR146KFX"]:
		
		encid=encode_id[name]
		model_paths = model_atac[model_atac[1]==name][2].values
		print(model_paths)
		
		if os.path.isfile(output_dir+"/"+encid+".json"):
			continue
		
		print(encid)

		args_json = {}
		
		success, args_json = main_fetch_preprocessing_files(encid, args_json, data_to_bam[name], name)
		if not success:
			print("ERR prep")
			continue

		success, args_json = main_fetch_bias_training_files(encid, args_json, model_paths, name)
		if not success:
			print("ERR bias prep")
			continue
		
		success, args_json = main_fetch_bias_model_files(encid, args_json, model_paths)
		if not success:
			print("ERR bias models")
			continue
			
		
		with open(output_dir+"/"+encid+".json", "w") as outfile:
			json.dump(args_json, outfile, indent=4)
	
	#print(args_json)

	
	
		
		
