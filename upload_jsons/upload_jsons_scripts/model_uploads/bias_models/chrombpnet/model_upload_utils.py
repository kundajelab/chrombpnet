import os
import json
import numpy as np


### utils for model uploads

def fetch_per_fold_models(odir, model_dir, encid, fold_num):
	input_paths = []
	log_paths = []
	log_paths_opt = []
	
	cmb = os.path.join(model_dir, "chrombpnet_model/chrombpnet_wo_bias.h5")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet_nobias.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		print(cmb)
		return None, None, None

	cmb = os.path.join(model_dir, "chrombpnet_model/chrombpnet.h5")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		print(cmb)
		return None, None, None
				
# 	checks_file = os.path.join(model_dir, "new_chrombpnet_model/check_passed.txt")
# 	if os.path.isfile(checks_file):
# 		cm_model = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet.h5")
# 		if os.path.isfile(cm_model):
# 			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".h5"))
# 		else:
# 			print(cm_model)
# 			return None, None, None
# 		
# 		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/chrombpnet.tar")
# 		if os.path.isfile(cm_model):
# 			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".tar"))
# 		else:
# 			print(cm_model)
# 			return None, None, None
# 
# 		
# 	else:
# 		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet_new.h5")
# 		if os.path.isfile(cm_model):
# 			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".h5"))
# 		else:
# 			print(cm_model)
# 			return None, None, None
# 		
# 		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar")
# 		if os.path.isfile(cm_model):
# 			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".tar"))
# 		else:
# 			print(cm_model)
# 			return None, None, None

						
	bm_model = os.path.join(model_dir, "chrombpnet_model/bias_model_scaled.h5")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias_scaled.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		print(cmb)
		return None, None, None

	cmb = os.path.join(model_dir, "new_model_formats_may_7_24_vf/chrombpnet.tar")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		print(cmb)

		return None, None, None

	cmb = os.path.join(model_dir, "new_model_formats_may_7_24_vf/chrombpnet_wo_bias.tar")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet_nobias.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		print(cmb)

		return None, None, None
		
			
	bm_model = os.path.join(model_dir, "new_model_formats_may_7_24_vf/bias_model_scaled.tar")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias_scaled.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None, None
		
	### fetch main logs
		
	modelling_log = os.path.join(model_dir, "chrombpnet_model/chrombpnet.args.json")
	if os.path.isfile(modelling_log):
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".args.json"))
	else:
		print(modelling_log)	
		
	modelling_log = os.path.join(model_dir, "chrombpnet_model/chrombpnet_data_params.tsv")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_data_params.tsv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(model_dir, "chrombpnet_model/chrombpnet_model_params.tsv")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_model_params.tsv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(model_dir, "chrombpnet_model/chrombpnet.params.json")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet.params.json"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(model_dir, "chrombpnet_model/chrombpnet.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".epoch_loss.csv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(model_dir, "chrombpnet_model/chrombpnet.log.batch")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".batch_loss.tsv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(model_dir, "chrombpnet_model/train_chrombpnet_model.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))
	else:
		print(modelling_log)
					
		
	return input_paths, log_paths, log_paths_opt
	
	
### utils for training and testing regions

def fetch_preprocessing_log_files(odir, encid, main_dir, name):
	# do bed file checks
	log_paths = []
	
	preprocessing_log = os.path.join(main_dir, name + "/data/"+name+"_preprocessing.log")
	if os.stat(preprocessing_log).st_size != 0:
			log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".stdout.txt"))

	try:
		preprocessing_log = os.path.join(main_dir, name + "/data/"+name.lower()+"_atac_fold_0.sh")
		if os.stat(preprocessing_log).st_size != 0:
			log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".script_v2.sh"))
	except:
		try:
			preprocessing_log = os.path.join(main_dir, name + "/data/"+name+"_DNASE_PE.sh")
			if os.stat(preprocessing_log).st_size != 0:
				log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".script_v2.sh"))
		except:
			preprocessing_log = os.path.join(main_dir, name + "/data/"+"h1_dnase_fold_0.sh")
			if os.stat(preprocessing_log).st_size != 0:
				log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".script_v2.sh"))
	
	preprocessing_log = os.path.join(main_dir, name + "/data/"+name+"_bias_pwm.png")
	if os.stat(preprocessing_log).st_size != 0:
			log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".bias_pwm.png"))

	return log_paths
	
def fetch_per_fold_training_data(odir,model_dir,encid, fold_num, main_dir, name):
	input_paths = []
	log_paths = []
	
	opath = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/splits_format/"
	filtered_regions_bed = os.path.join(opath + "/fold_"+str(fold_num)+".json")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"cv_params.fold_"+str(fold_num)+".json"))

	if fold_num==0:
		filtered_regions_bed = os.path.join(main_dir, name+"/negatives_data/negatives_with_summit.bed.gz")
		if os.path.isfile(filtered_regions_bed):
			input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))
	else:
		filtered_regions_bed = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/negatives_with_summit.bed.gz")
		if os.path.isfile(filtered_regions_bed):
			input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))


	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_may_7_2024/peaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_may_7_2024/peaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_may_7_2024/peaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_may_7_2024/nonpeaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_may_7_2024/nonpeaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_may_7_2024/nonpeaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))
			
	# preprocessing logs to include

	if fold_num==0:
		#negatives_log = os.path.join(temp_dir, name+"/negatives_data/make_background_regions.log")
		negatives_log = os.path.join(main_dir, name+"/negatives_data/make_background_regions.log")

		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))
	else:
		negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/make_background_regions.log")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))


	if fold_num==0:
#		negatives_log = os.path.join(temp_dir, "negatives_data/negatives_compared_with_foreground.png")
		negatives_log = os.path.join(main_dir, name+"/negatives_data/negatives_compared_with_foreground.png")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))
	else:
		negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/negatives_compared_with_foreground.png")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))
	
	return input_paths, log_paths
