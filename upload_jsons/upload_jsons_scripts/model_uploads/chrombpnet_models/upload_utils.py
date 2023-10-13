import os
import json
import numpy as np

### utils for preprocessing

def fetch_input_bam_ids(odir,encid):
	log_path = os.path.join(odir, encid + "/preprocessing/preprocess_"+encid+".log")
	logd = open(log_path).readlines()
	set_cflag=False
	set_bflag=False
	
	bams_ids = []
	
	for line in logd:
	
		if set_cflag:
			words = line.strip().split()
			if words[1] == "cp":
				if words[2].split("/")[-1].endswith("bam"):
					bam_enc = words[2].split("/")[-1].replace(".bam","")
					bams_ids.append(bam_enc)
					return bams_ids
				else:
					print(encid,"error")
					return
			else:
				print(encid,"error")
				return
		
		if set_bflag:
			words = line.strip().split()
			if words[1] == "samtools" and words[2] == "merge":
				encids = words[6:]
				for encid in encids:
					if encid.split("/")[-1].endswith(".bam"):
						bam_enc = encid.split("/")[-1].replace(".bam","")
						bams_ids.append(bam_enc)
					else:
						print(encid,"error")
						return
				return bams_ids
			else:
				print(encid,"error")
				return
				
		if "Only one source bam file found. Copying over as merged file." in line:
			set_cflag=True
		if "Merging bam files" in line:
			set_bflag=True	
		
### utils for training and testing regions

def fetch_preprocessing_log_files(odir, encid):
	# do bed file checks
	log_paths = []
	
	# preprocessing (6 files)
	preprocessing_log = os.path.join(odir, encid + "/preprocessing/preprocessing.log.e")
	if os.stat(preprocessing_log).st_size != 0:
		log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".stderr.txt"))

	preprocessing_log = os.path.join(odir, encid + "/preprocessing/preprocessing.log.o")
	if os.stat(preprocessing_log).st_size != 0:
		log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".stdout.txt"))

	preprocessing_log = os.path.join(odir, encid + "/preprocessing/"+encid+".log")
	if os.stat(preprocessing_log).st_size != 0:
		log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".stdout_v1.txt"))

	preprocessing_log = os.path.join(odir, encid + "/preprocessing/preprocess_"+encid+".log")
	if os.stat(preprocessing_log).st_size != 0:
		log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".stdout_v2.txt"))

	preprocessing_log = os.path.join(odir, encid + "/preprocessing/params_file.json")
	if os.stat(preprocessing_log).st_size != 0:
		log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".params_file.json"))

	preprocessing_log = os.path.join(odir, encid + "/preprocessing/bigWigs/"+encid+".png")
	if os.stat(preprocessing_log).st_size != 0:
			log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".bias_pwm.png"))

	# peak_logs (2 files)
	negatives_log = os.path.join(odir, encid + "/peak_calling/log.e")
	if os.path.isfile(negatives_log):
		log_paths.append((negatives_log,"logfile.peak_calling."+encid+".stdout_v1.txt"))

	negatives_log = os.path.join(odir, encid + "/peak_calling/log.o")
	if os.path.isfile(negatives_log):
		log_paths.append((negatives_log,"logfile.peak_calling."+encid+".stdout_v2.txt"))

	# negative logs	(4 files)
	negatives_log = os.path.join(odir, encid + "/negatives_data/make_background_regions.log")
	if os.path.isfile(negatives_log):
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching."+encid+".stdout_v1.txt"))

	negatives_log = os.path.join(odir, encid + "/negatives_data/"+encid+".log")
	if os.path.isfile(negatives_log):
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching."+encid+".stdout_v1.txt"))


	negatives_log = os.path.join(odir, encid + "/negatives_data/gc_matching.log.o")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.gc_matching."+encid+".stdout.txt"))
		
	negatives_log = os.path.join(odir, encid + "/negatives_data/gc_matching.log.e")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.gc_matching."+encid+".stderr.txt"))


	negatives_log = os.path.join(odir, encid + "/negatives_data/negatives_compared_with_foreground.png")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.gc_matching."+encid+".stdout.png"))

	# all test log files include once (+ 2 files)
	try:
		negatives_log = os.path.join(odir, encid + "/negatives_data/test/test.gc_matching.log.o")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.test.gc_matching."+encid+".stdout.txt"))

		negatives_log = os.path.join(odir, encid + "/negatives_data/test/test.gc_matching.log.e")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.test.gc_matching."+encid+".stderr.txt"))

	except:
		negatives_log = os.path.join(odir, encid + "/negatives_data/test/test..gc_matching.log.o")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.test.gc_matching."+encid+".stdout.txt"))

		negatives_log = os.path.join(odir, encid + "/negatives_data/test/test..gc_matching.log.e")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.test.gc_matching."+encid+".stderr.txt"))
	

	return log_paths
	
def fetch_per_fold_training_data(odir,model_dir,encid, fold_num):
	input_paths = []
	log_paths = []
	
	opath = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/splits_format/"
	filtered_regions_bed = os.path.join(opath + "/fold_"+str(fold_num)+".json")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"cv_params.fold_"+str(fold_num)+".json"))

	filtered_regions_bed = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/peaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/peaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/peaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/nonpeaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/nonpeaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/nonpeaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))
			
	# preprocessing logs to include

	negatives_log = os.path.join(odir, encid + "/negatives_data/test/fold_"+str(fold_num)+"."+encid+"_test.log")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.test.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	negatives_log = os.path.join(odir, encid + "/negatives_data/test/test.fold_"+str(fold_num)+".negatives_compared_with_foreground.png")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.test.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))

	negatives_log = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/pconv.log.o")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.regions.fold_"+str(fold_num)+"."+encid+".formatting.stdout.txt"))

	negatives_log = os.path.join(odir, encid + "/" + model_dir + "/train_test_regions/pconv.log.e")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.regions.fold_"+str(fold_num)+"."+encid+".formatting.stderr.txt"))
		
	return input_paths, log_paths

### utils for model uploads

def fetch_per_fold_models(odir, model_dir, encid, fold_num):
	input_paths = []
	log_paths = []
	log_paths_opt = []
	
	cmb = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet_wo_bias.h5")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet_nobias.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		return None, None, None
		
	checks_file = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/check_passed.txt")
	if os.path.isfile(checks_file):
		cm_model = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet.h5")
		if os.path.isfile(cm_model):
			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".h5"))
		else:
			print(cm_model)
			return None, None, None
		
		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/chrombpnet.tar")
		if os.path.isfile(cm_model):
			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".tar"))
		else:
			print(cm_model)
			return None, None, None

		modelling_log = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/conv.chrombpnet.log.e")
		if os.stat(modelling_log).st_size != 0:
			log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_formatting.stderr.txt"))
		else:
			print(modelling_log)
		
		modelling_log = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/conv.chrombpnet.log.o")
		if os.stat(modelling_log).st_size != 0:
			log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_formatting.stdout.txt"))
		else:
			print(modelling_log)

		
	else:
		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet_new.h5")
		if os.path.isfile(cm_model):
			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".h5"))
		else:
			print(cm_model)
			return None, None, None
		
		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar")
		if os.path.isfile(cm_model):
			input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".tar"))
		else:
			print(cm_model)
			return None, None, None

		modelling_log = os.path.join(odir, encid + "/fix_h5_logs/fix.h5.log.e")
		if os.stat(modelling_log).st_size != 0:
			log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_fix_formatting.stderr.txt"))
		else:
			print(modelling_log)
		
		modelling_log = os.path.join(odir, encid + "/fix_h5_logs/fix.h5.log.o")
		if os.stat(modelling_log).st_size != 0:
			log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_fix_formatting.stdout.txt"))
		else:
			print(modelling_log)

						
	bm_model = os.path.join(odir, encid + "/" + model_dir + "/bias_model_scaled.h5")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias_scaled.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		return None, None, None

	cmb = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/chrombpnet_wo_bias.tar")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet_nobias.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None, None
		
			
	bm_model = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/bias_model_scaled.tar")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias_scaled.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None, None
		
	### fetch main logs
		
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet.args.json")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".args.json"))
	else:
		print(modelling_log)	
		
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet_data_params.tsv")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_data_params.tsv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet_model_params.tsv")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_model_params.tsv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet.params.json")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet.params.json"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".epoch_loss.csv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/chrombpnet.log.batch")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".batch_loss.tsv"))
	else:
		print(modelling_log)
		
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/train_chrombpnet_model.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))
	else:
		print(modelling_log)
					
	#### fetch model training log files ########
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/modelling.log.e")
	if os.path.isfile(modelling_log):
		if os.stat(modelling_log).st_size != 0:
			log_paths_opt.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".stderr.txt"))
		else:
			print(modelling_log)
	else:
		print(modelling_log)
			
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/modelling.log.o")
	if os.path.isfile(modelling_log):
		if os.stat(modelling_log).st_size != 0:
			log_paths_opt.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".stdout.txt"))
		else:
			print(modelling_log)
	else:
		print(modelling_log)
		
	#### fetch model conversion log files ########
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/conv.bias_model_scaled.log.e")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".bias_formatting.stderr.txt"))
	else:
		print(modelling_log)
			
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/conv.bias_model_scaled.log.o")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".bias_formatting.stdout.txt"))
	else:
		print(modelling_log)
				
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/conv.chrombpnet_wo_bias.log.o")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_no_bias_formatting.stderr.txt"))
	else:
		print(modelling_log)
			
	modelling_log = os.path.join(odir, encid + "/" + model_dir + "/new_model_formats/conv.chrombpnet_wo_bias.log.o")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".chrombpnet_no_bias_formatting.stdout.txt"))
	else:
		print(modelling_log)
		
	return input_paths, log_paths, log_paths_opt




			
			
	
