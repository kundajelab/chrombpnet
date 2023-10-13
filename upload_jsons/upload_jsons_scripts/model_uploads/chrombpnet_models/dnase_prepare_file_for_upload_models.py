import os
import upload_utils
import json

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"
bw_odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/full_deepshaps/bigwigs/DNASE/"
#output_dir =  "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/DNASE/stage1/jul_26_2023/"
output_dir="dnase_production_uploads/"

tissue_encids = open("data/tissue_passed.txt").readlines()
tissue_encids = [line.strip() for line in tissue_encids]

primary_encids = open("data/primary_passed.txt").readlines()
primary_encids = [line.strip() for line in primary_encids]

celline_encids = open("data/cellline_passed.txt").readlines()
celline_encids = [line.strip() for line in celline_encids]

invitro_encids = open("data/invitro_passed.txt").readlines()
invitro_encids = [line.strip() for line in invitro_encids]

primary_models_path = ["chrombppnet_model_encsr283tme_bias", "chrombppnet_model_encsr283tme_bias_fold_1", "chrombppnet_model_encsr283tme_bias_fold_2", "chrombppnet_model_encsr283tme_bias_fold_3", "chrombppnet_model_encsr283tme_bias_fold_4"]
celline_models_path = ["chrombpnet_model_feb15_fold_0", "chrombpnet_model_feb15_fold_1", "chrombpnet_model_feb15_fold_2", "chrombpnet_model_feb15_fold_3", "chrombpnet_model_feb15_fold_4"]
tissue_models_path = ["chrombpnet_model_encsr880cub_bias","chrombppnet_model_encsr880cub_bias_fold_1","chrombppnet_model_encsr880cub_bias_fold_2","chrombppnet_model_encsr880cub_bias_fold_3","chrombppnet_model_encsr880cub_bias_fold_4"]
invitro_models_path = ["chrombpnet_model_encsr146kfx_bias", "chrombpnet_model_encsr146kfx_bias_fold_1", "chrombpnet_model_encsr146kfx_bias_fold_2", "chrombpnet_model_encsr146kfx_bias_fold_3", "chrombpnet_model_encsr146kfx_bias_fold_4"]

encids = tissue_encids + primary_encids + celline_encids + invitro_encids

def main_fetch_preprocessing_files(encid, args_json, bias_encid):

	success_flag = False
	args_json["upload bias"] = False
	args_json["bias model encid"] = bias_encid

	# find the bams input
	preprocessing_path = os.path.join(odir, encid + "/preprocessing/bigWigs/"+encid+".bigWig")
	if os.path.isfile(preprocessing_path):
		bam_ids = upload_utils.fetch_input_bam_ids(odir,encid)
		
		if bam_ids == None:
			success = False
			return  success_flag, args_json
			
		args_json["experiment"] = encid
		args_json["bam files"] = bam_ids
		args_json["assay"] = "DNase-seq"
		args_json["observed signal profile bigWig"] = preprocessing_path
		success = True
	else:
		success = False
	
	return	success, args_json

def main_fetch_model_files(encid, args_json):
	success = False
	args_json["models tar"] = {}
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/models.README"
	assert(os.path.isfile(readme_file))
	args_json["models tar"]["file.paths"] = [(readme_file, "README.md")]
	#args_json["models tar"]["logs.models."+encid] = {"file.paths": None}

	for i in range(5):
		data_paths, log_paths, log_paths_opt = upload_utils.fetch_per_fold_models(odir,models_path[i], encid, i)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["models tar"]["fold_"+str(i)] = {}
		args_json["models tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["models tar"]["fold_"+str(i)]["logs.models.fold_"+str(i)+"."+encid] = {"file.paths": log_paths+log_paths_opt}
		# 9 log file expected per model
		assert(len(data_paths) == 6)
		assert(len(log_paths) == 13)

	success=True
	return success, args_json

def main_fetch_training_files(encid, args_json):
	success = False
	
	# find the training test regions
	args_json["training and test regions tar"] = {}	
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/training_test_regions.README"
	assert(os.path.isfile(readme_file))
	args_json["training and test regions tar"]["file.paths"] = [(readme_file, "README.md")]

	input_peaks = os.path.join(odir, encid + "/preprocessing/downloads/peaks.bed.gz")
	if os.path.isfile(input_peaks):
		args_json["training and test regions tar"]["file.paths"].append((input_peaks,"peaks.all_input_regions."+encid+".bed.gz"))		
	else:
		success = False
		return success, args_json

	input_nonpeaks_gz = os.path.join(odir, encid + "/negatives_data/negatives_with_summit.bed.gz")
	input_nonpeaks = os.path.join(odir, encid + "/negatives_data/negatives_with_summit.bed")
	if not os.path.isfile(input_nonpeaks_gz):
		if os.path.isfile(input_nonpeaks):
			import pandas as pd
			#os.system("gzip "+input_nonpeaks)
			nonpeaks_data = pd.read_csv(input_nonpeaks, sep="\t", header=None)
			nonpeaks_data.to_csv(input_nonpeaks+".gz", sep="\t", header=False, index=False, compression="gzip")
		#os.system("rm "+input_nonpeaks)

	input_nonpeaks = os.path.join(odir, encid + "/negatives_data/negatives_with_summit.bed.gz")

	if os.path.isfile(input_nonpeaks):
		args_json["training and test regions tar"]["file.paths"].append((input_nonpeaks,"nonpeaks.all_input_regions."+encid+".bed.gz"))
	else:
		success = False
		return success, args_json

	log_paths = upload_utils.fetch_preprocessing_log_files(odir,encid)
	args_json["training and test regions tar"]["logs.training_test_regions."+encid] = {"file.paths": log_paths}
	#print(len(log_paths))
	#print(log_paths)
	assert(len(log_paths) == 14)		

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_training_data(odir,models_path[i], encid, i)

		args_json["training and test regions tar"]["fold_"+str(i)] = {}
		args_json["training and test regions tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["training and test regions tar"]["fold_"+str(i)]["logs.training_test_regions.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		#print(data_paths)
		assert(len(data_paths) == 7)
		
		assert(len(log_paths) == 4)	

		if len(data_paths) != 7:
			success = False
			return success, args_json
	
	success = True
	return success, args_json

		
if __name__ == "__main__":

	ignore_list = []
	# missing args file
	ignore_list += ["ENCSR720TCN", "ENCSR241BNZ", "ENCSR516JCM", "ENCSR642DZF", "ENCSR787ERP"]
	# missing train_test
	ignore_list += ['ENCSR000EOT']
	# chrombpnet models
	ignore_list += ["ENCSR000EMT", "ENCSR149XIL", "ENCSR477RTP", "ENCSR000EOT", "ENCSR000EMU"]
	
	for encid in encids:
		if encid  in ignore_list:
			continue
		
		if encid in primary_encids:
			models_path = primary_models_path
			bias_encid="ENCSR283TME"
			#print("primary")
		elif encid in tissue_encids:
			models_path = tissue_models_path
			bias_encid="ENCSR880CUB"
			#print("tissue")
		elif encid in invitro_encids:
			models_path = invitro_models_path
			bias_encid="ENCSR146KFX"
			#print("invitro")
		elif encid in celline_encids:
			models_path = celline_models_path
			bias_encid="ENCSR149XIL"
			#print("celline")
		else:
			print(encid)
			print("type not found")
			continue
				
		if os.path.isfile(output_dir+"/"+encid+".json"):
			continue
		
		print(encid)
		args_json = {}
		
		
		success, args_json = main_fetch_preprocessing_files(encid, args_json, bias_encid)
		if not success:
			print(encid)
			print("exit preprocessing")
			continue
		
		success, args_json = main_fetch_model_files(encid, args_json)
		if not success:
			print(encid)
			print("exit models")
			continue

		success, args_json = main_fetch_training_files(encid, args_json)
		if not success:
			print(encid)
			print("exit train test regions")
			continue

		
		with open(output_dir+"/"+encid+".json", "w") as outfile:
			json.dump(args_json, outfile, indent=4)
	
	#print(args_json)

	
	
		
		
