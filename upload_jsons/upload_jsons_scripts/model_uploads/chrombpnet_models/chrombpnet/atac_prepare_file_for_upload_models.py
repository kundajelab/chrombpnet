import os
import upload_utils
import json

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/"
bw_odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/full_deepshaps/bigwigs/ATAC/"
#output_dir =  "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/ATAC/stage1/jul_17_2023/"
models_path = ["chrombpnet_model_feb15", "chrombpnet_model_feb15_fold_1", "chrombpnet_model_feb15_fold_2", "chrombpnet_model_feb15_fold_3", "chrombpnet_model_feb15_fold_4"]
output_dir = "atac_production_uploads/"
#encids = os.listdir(odir)
encids = open("data/atac_passed.txt").readlines()
encids = [line.strip() for line in encids]

encode_id = {"K562": "ENCSR868FGK",
"GM12878": "ENCSR637XSC",
"HEPG2": "ENCSR291GJU",
"IMR90": "ENCSR200OML",
"H1ESC": "GSE267154"}     

def main_fetch_preprocessing_files(encid, args_json):

	success_flag = False
	args_json["upload bias"] = False
	args_json["bias model encid"] = encid

	# find the bams input
	preprocessing_path = os.path.join(odir, encid + "/preprocessing/bigWigs/"+encid+".bigWig")
	if os.path.isfile(preprocessing_path):
		bam_ids = upload_utils.fetch_input_bam_ids(odir,encid)
		
		if bam_ids == None:
			success = False
			return  success_flag, args_json
			
		args_json["experiment"] = encid
		args_json["bam files"] = bam_ids
		args_json["assay"] = "ATAC-seq"
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
	args_json["models tar"]["logs.models."+encid] = {"file.paths": None}

	for i in range(5):
		data_paths, log_paths, log_paths_opt = upload_utils.fetch_per_fold_models(odir,models_path[i], encid, i)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["models tar"]["fold_"+str(i)] = {}
		args_json["models tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["models tar"]["fold_"+str(i)]["logs.models.fold_"+str(i)+"."+encid] = {"file.paths": log_paths+log_paths_opt}
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
		
	input_nonpeaks = os.path.join(odir, encid + "/negatives_data/negatives_with_summit.bed")
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
	assert(len(log_paths) == 12)		

	for i in range(5):
		data_paths, log_paths = upload_utils.fetch_per_fold_training_data(odir,models_path[i], encid, i)

		args_json["training and test regions tar"]["fold_"+str(i)] = {}
		args_json["training and test regions tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["training and test regions tar"]["fold_"+str(i)]["logs.training_test_regions.fold_"+str(i)+"."+encid] = {"file.paths": log_paths}
		assert(len(data_paths) == 7)
		assert(len(log_paths) == 4)		

		if len(data_paths) != 7:
			success = False
			return success, args_json
	
	success = True
	return success, args_json
	
		
if __name__ == "__main__":


	for name in ["K562", "GM12878", "HEPG2", "IMR90", "H1ESC"]:

		
		encid=encode_id[name]
		if os.path.isfile(output_dir+"/"+encid+".json"):
			continue
		
		print(encid)

		args_json = {}
		
		success, args_json = main_fetch_preprocessing_files(encid, args_json)
		if not success:
			print("fail prep")
			continue
		
		success, args_json = main_fetch_model_files(encid, args_json)
		if not success:
			print("fail model")
			continue

		success, args_json = main_fetch_training_files(encid, args_json)
		if not success:
			print("fail train prep")
			continue

		
		with open(output_dir+"/"+encid+".json", "w") as outfile:
			json.dump(args_json, outfile, indent=4)
	
	#print(args_json)

	
	
		
		
