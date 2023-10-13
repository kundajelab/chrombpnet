import os
import upload_utils
import json
import pandas as pd

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/"
output_dir =  "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/ATAC/stage1/jul_17_2023/"
main_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"
output_dir = "temp/"

encids = os.listdir(odir)
#encids = open("../chromatin_atlas_atac/test_encid.txt").readlines()
#encids = [line.strip() for line in encids]

model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",", header=None)
encode_id = {"K562": "ENCSR868FGK",
        "GM12878": "ENCSR637XSC",
        "HEPG2": "ENCSR291GJU",
        "IMR90": "ENCSR200OML"}
        
data_to_bam = {"K562": ["ENCFF077FBI", "ENCFF128WZG", "ENCFF534DCE"],
        "GM12878": ["ENCFF440GRZ", "ENCFF962FMH", "ENCFF981FXV"],
        "HEPG2": ["ENCFF926KFU", "ENCFF624SON", "ENCFF990VCP"],
        "IMR90": ["ENCFF848XMR", "ENCFF715NAV"]}
        
def main_fetch_preprocessing_files(encid, args_json, bam_ids, name):
	# define bam_ids, name
	
	success_flag = False
	
	args_json["upload bias"] = True
	args_json["bias model encid"] = encid 

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

def main_fetch_bias_model_files(encid, args_json, models_path):
	success = False
	args_json["bias models tar"] = {}
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/models.README"
	assert(os.path.isfile(readme_file))
	args_json["bias models tar"]["file.paths"] = [(readme_file, "README.md")]
	args_json["bias models tar"]["logs.bias.models."+encid] = {"file.paths": None}

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_bias_models(odir, models_path[i], encid, i)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["bias models tar"]["fold_"+str(i)] = {}
		args_json["bias models tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["bias models tar"]["fold_"+str(i)]["logs.bias.models.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		# 9 log file expected per model
		assert(len(log_paths) == 2)	
		assert(len(data_paths) == 2)		
	success=True
	return success, args_json
	
def main_fetch_model_files(encid, args_json, models_path):
	success = False
	args_json["models tar"] = {}
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/models.README"
	assert(os.path.isfile(readme_file))
	args_json["models tar"]["file.paths"] = [(readme_file, "README.md")]
	args_json["models tar"]["logs.models."+encid] = {"file.paths": None}

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_models(odir, models_path[i], encid, i)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["models tar"]["fold_"+str(i)] = {}
		args_json["models tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["models tar"]["fold_"+str(i)]["logs.models.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		# 9 log file expected per model
		assert(len(log_paths) == 2)		
	success=True
	return success, args_json


def main_fetch_bias_training_files(encid, args_json, models_path, name):
	success = False
	
	# find the training test regions
	args_json["bias training and test regions tar"] = {}	
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/training_test_regions.README"
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
	assert(len(log_paths) == 2)		

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_training_data_bias(odir, models_path[i], encid, i, main_dir, name)

		args_json["bias training and test regions tar"]["fold_"+str(i)] = {}
		args_json["bias training and test regions tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["bias training and test regions tar"]["fold_"+str(i)]["logs.bias.training_test_regions.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		assert(len(data_paths) == 6)
		assert(len(log_paths) == 3)	

		#if len(data_paths) != 3:
		#	success = False
		#	return success, args_json
	
	success = True
	return success, args_json
	
def main_fetch_training_files(encid, args_json, models_path, name):
	success = False
	
	# find the training test regions
	args_json["training and test regions tar"] = {}	
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/training_test_regions.README"
	assert(os.path.isfile(readme_file))
	args_json["training and test regions tar"]["file.paths"] = [(readme_file, "README.md")]

	input_peaks = os.path.join(main_dir, name + "/data/peaks_no_blacklist.bed.gz")
	if os.path.isfile(input_peaks):
		args_json["training and test regions tar"]["file.paths"].append((input_peaks,"peaks.all_input_regions."+encid+".bed.gz"))		
	else:
		success = False
		return success, args_json
		
	log_paths = upload_utils.fetch_preprocessing_log_files(odir, encid, main_dir, name)
	args_json["training and test regions tar"]["logs.training_test_regions."+encid] = {"file.paths": log_paths}
	assert(len(log_paths) == 4)		

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_training_data(odir, models_path[i], encid, i, main_dir, name)

		args_json["training and test regions tar"]["fold_"+str(i)] = {}
		args_json["training and test regions tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["training and test regions tar"]["fold_"+str(i)]["logs.training_test_regions.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		assert(len(data_paths) == 8)
		assert(len(log_paths) == 3)	

		#if len(data_paths) != 3:
		#	success = False
		#	return success, args_json
	
	success = True
	return success, args_json

def main_fetch_perf_metric_files(encid, args_json, models_path, name):

	args_json["model performance metrics tar"] = {}
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/READMES/performance_metrics.README"
	assert(os.path.isfile(readme_file))
	args_json["model performance metrics tar"]["file.paths"] = [(readme_file, "README.md")]

	temp_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/ATAC_PE/"
	input_png = os.path.join(temp_dir, name + "/data/"+name+"_bias_pwm.png")
	if os.path.isfile(input_png):
		args_json["model performance metrics tar"]["file.paths"].append((input_png,"bias_motif.shift_qc.obs_signal_profile."+encid+".png"))		
	else:
		success = False
		return success, args_json

	#temp_dir="/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/"
	#input_sam_qc = os.path.join(temp_dir, name + "/" + name+".samstats.qc")
	#if os.path.isfile(input_sam_qc):
	#	args_json["model performance metrics tar"]["file.paths"].append((input_sam_qc,"samstats.obs_signal_profile."+encid+".txt"))		
	#else:
	#	success = False
	#	return success, args_json

	log_paths = upload_utils.fetch_model_perf_logs(odir, encid, name)
	args_json["model performance metrics tar"]["logs.performance_metrics."+encid] = {"file.paths": log_paths}
	#print(log_paths)
	assert(len(log_paths) == 1)
	
	
	for i in range(5):	
		data_paths, log_paths, footprint_paths = upload_utils.fetch_per_fold_model_metrics(odir,models_path[i], encid, i, "ATAC")

		args_json["model performance metrics tar"]["fold_"+str(i)] = {}
		args_json["model performance metrics tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["model performance metrics tar"]["fold_"+str(i)]["logs.performance_metrics.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		args_json["model performance metrics tar"]["fold_"+str(i)]["footprint_plot.bias_motif.fold_"+str(i)+"."+encid] = {"file.paths": footprint_paths}

		print(data_paths)
		print(len(data_paths))
		assert(len(data_paths) >= 12)
		#assert(len(log_paths) == 2)
		#assert(len(footprint_paths) == 6)
		
		if len(data_paths) < 12:
			success = False
			return success, args_json

		#if len(footprint_paths) != 6:
		#	success = False
		#	return success, args_json
					
	success, median_json_chrombpnet =  upload_utils.fetch_chrombpnet_median_performance(odir, encid, args_json, models_path, main_dir, name)
	if success == False:
		return success, args_json
	args_json["bpnet quality metrics json"] = median_json_chrombpnet
	args_json["model performance metrics tar"]["file.paths"].append((median_json_chrombpnet, "performance_metrics.chrombpnet.fold_median."+encid+".json"))
	success, median_json_chrombpnet_nobias =  upload_utils.fetch_chrombpnet_nobias_median_performance(odir, encid, args_json, models_path, main_dir, name)
	if success == False:
		return success, args_json
	args_json["model performance metrics tar"]["file.paths"].append((median_json_chrombpnet_nobias, "performance_metrics.chrombpnet_nobias.fold_median."+encid+".json"))

	
	success=True
	return success, args_json
	
		
if __name__ == "__main__":

	for name in ["K562", "GM12878", "HEPG2"]:
		
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
		
		#success, args_json = main_fetch_model_files(encid, args_json, model_paths)
		#if not success:
		#	print("ERR models")
		#	continue

		success, args_json = main_fetch_bias_model_files(encid, args_json, model_paths)
		if not success:
			print("ERR bias models")
			continue
			
		#success, args_json = main_fetch_training_files(encid, args_json, model_paths, name)
		#if not success:
		#	print("ERR training data")
		#	continue
			
		#success, args_json = main_fetch_perf_metric_files(encid, args_json, model_paths, name)
		#if not success:
		#	print("ERR model perf")
		#	continue
		
		with open(output_dir+"/"+encid+".json", "w") as outfile:
			json.dump(args_json, outfile, indent=4)
	
	#print(args_json)

	
	
		
		
