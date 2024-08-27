import os
import atac_bias_upload_utils as upload_utils
import json
import pandas as pd
import model_upload_utils

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/"
#output_dir =  "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/ATAC/stage1/jul_17_2023/"
main_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"
output_dir = "atac_production_uploads/"

encids = os.listdir(odir)
#encids = open("../chromatin_atlas_atac/test_encid.txt").readlines()
#encids = [line.strip() for line in encids]

model_atac = pd.read_csv("atac_bias_model_chrombpnet.csv",sep=",", header=None)
encode_id = {"K562": "ENCSR868FGK",
"GM12878": "ENCSR637XSC",
"HEPG2": "ENCSR291GJU",
"IMR90": "ENCSR200OML",
"H1ESC": "GSE267154"}     
  
data_to_bam = {"K562": ["ENCFF077FBI", "ENCFF128WZG", "ENCFF534DCE"],
"GM12878": ["ENCFF440GRZ", "ENCFF981FXV", "ENCFF962FMH"],
"HEPG2": ["ENCFF624SON", "ENCFF926KFU", "ENCFF990VCP"],
"IMR90": ["ENCFF848XMR", "ENCFF715NAV"],
"H1ESC": ["GSM8260976", "GSM8260977"]
} 

def main_fetch_training_files(encid, args_json, model_paths, name):
	success = False
	
	# find the training test regions
	args_json["training and test regions tar"] = {}	
	readme_file = "READMEs/training.README"
	assert(os.path.isfile(readme_file))
	args_json["training and test regions tar"]["file.paths"] = [(readme_file, "README.md")]

	input_peaks = os.path.join(main_dir, name + "/data/peaks_no_blacklist.bed.gz")
	if os.path.isfile(input_peaks):
		args_json["training and test regions tar"]["file.paths"].append((input_peaks,"peaks.all_input_regions."+encid+".bed.gz"))		
	else:
		success = False
		return success, args_json
		
	log_paths = model_upload_utils.fetch_preprocessing_log_files(odir,encid,main_dir, name)
	args_json["training and test regions tar"]["logs.training_test_regions."+encid] = {"file.paths": log_paths}
	assert(len(log_paths) == 3)		

	for i in range(5):
		data_paths, log_paths = model_upload_utils.fetch_per_fold_training_data(odir,model_paths[i], encid, i, main_dir, name)

		args_json["training and test regions tar"]["fold_"+str(i)] = {}
		args_json["training and test regions tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["training and test regions tar"]["fold_"+str(i)]["logs.training_test_regions.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		print(len(data_paths))
		assert(len(data_paths) == 8)
		assert(len(log_paths) == 2)		

		if len(data_paths) != 8:
			success = False
			return success, args_json
	
	success = True
	return success, args_json
	
def main_fetch_preprocessing_files_for_k562(encid, args_json, bam_ids, name):
	# define bam_ids, name
	
	success_flag = False
	
	args_json["upload bias"] = False
	#args_json["bias model encid"] = encid 

	# find the bams input
	preprocessing_path = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/data/"+name+"_unstranded.bw"

	if os.path.isfile(preprocessing_path):	
		args_json["experiment"] = encid
		args_json["bam files"] = bam_ids
		args_json["assay"] = "ATAC-seq"
		success = True
	else:
		success = False
	
	return	success, args_json
	
def main_fetch_preprocessing_files(encid, args_json, bam_ids, name):
	# define bam_ids, name
	
	success_flag = False
	
	args_json["upload bias"] = True
	#args_json["bias model encid"] = encid 

	# find the bams input
	preprocessing_path = "/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"+name+"/data/"+name+"_unstranded.bw"

	if os.path.isfile(preprocessing_path):	
		args_json["experiment"] = encid
		args_json["bam files"] = bam_ids
		args_json["assay"] = "ATAC-seq"
		args_json["observed signal profile bigWig"] = preprocessing_path
		success = True
	else:
		success = False
	
	return	success, args_json

def main_fetch_model_files(encid, args_json, model_paths, name):
	success = False
	args_json["models tar"] = {}
	readme_file = "READMEs/models.README"
	assert(os.path.isfile(readme_file))
	args_json["models tar"]["file.paths"] = [(readme_file, "README.md")]
	args_json["models tar"]["logs.models."+encid] = {"file.paths": None}

	for i in range(5):
		data_paths, log_paths, log_paths_opt = model_upload_utils.fetch_per_fold_models(odir,model_paths[i], encid, i)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["models tar"]["fold_"+str(i)] = {}
		args_json["models tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["models tar"]["fold_"+str(i)]["logs.models.fold_"+str(i)+"."+encid] = {"file.paths": log_paths+log_paths_opt}
		assert(len(data_paths) == 6)
		print(len(log_paths))
		assert(len(log_paths) >= 6)
				
	success=True
	return success, args_json
	
def main_fetch_bias_model_files(encid, args_json, models_path):
	success = False
	args_json["bias models tar"] = {}
	readme_file = "READMEs/bias.models.README" 
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
		# 9 log file expected per model
		print(len(log_paths))
		assert(len(log_paths) >= 2)	
		assert(len(data_paths) == 2)		
	success=True
	return success, args_json
	
def main_fetch_bias_training_files(encid, args_json, models_path, name):
	success = False
	
	# find the training test regions
	args_json["bias training and test regions tar"] = {}	
	readme_file = "READMEs/bias.training.README"
	assert(os.path.isfile(readme_file))
	args_json["bias training and test regions tar"]["file.paths"] = [(readme_file, "README.md")]

	input_peaks = os.path.join(main_dir, name + "/data/peaks_no_blacklist.bed.gz")
	print(input_peaks)
	if os.path.isfile(input_peaks):
		args_json["bias training and test regions tar"]["file.paths"].append((input_peaks,"peaks.all_input_regions."+encid+".bed.gz"))		
	else:
		success = False
		return success, args_json
		
	log_paths = upload_utils.bias_fetch_preprocessing_log_files(odir, encid, main_dir, name)
	args_json["bias training and test regions tar"]["logs.bias.training_test_regions."+encid] = {"file.paths": log_paths}
	assert(len(log_paths) == 3)

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_training_data_bias(odir, models_path[i], encid, i, main_dir, name)

		args_json["bias training and test regions tar"]["fold_"+str(i)] = {}
		args_json["bias training and test regions tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["bias training and test regions tar"]["fold_"+str(i)]["logs.bias.training_test_regions.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		#print(len(data_paths))
		#print(data_paths)
		assert(len(data_paths) == 5)
		assert(len(log_paths) == 2)	

		#if len(data_paths) != 3:
		#	success = False
		#	return success, args_json
	
	success = True
	return success, args_json
	
	
		
if __name__ == "__main__":

	for name in ["K562", "GM12878", "HEPG2", "IMR90", "H1ESC"]:
		
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
		
		if name == "K562":
			with open(output_dir+"/"+encid+"_bias.json", "w") as outfile:
				json.dump(args_json, outfile, indent=4)
			
			args_json = {}
			main_fetch_preprocessing_files_for_k562(encid, args_json, data_to_bam[name], name)
			if not success:
				print("ERR prep")
				continue
				
		success, args_json = main_fetch_model_files(encid, args_json, model_paths, name)
		if not success:
			print("fail model")
			continue
		
		success, args_json = main_fetch_training_files(encid, args_json, model_paths, name)
		if not success:
			print("fail train prep")
			continue

	
		with open(output_dir+"/"+encid+".json", "w") as outfile:
			json.dump(args_json, outfile, indent=4)
	
	#print(args_json)

	
	
		
		
