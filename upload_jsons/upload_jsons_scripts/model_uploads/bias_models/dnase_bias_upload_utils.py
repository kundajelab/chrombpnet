import os
import json
import numpy as np

### utils for preprocessing

		
### utils for training and testing regions


def bias_fetch_preprocessing_log_files_set_1(odir, encid, main_dir, name):
	# do bed file checks
	log_paths = []
	
	
	# preprocessing log files
	preprocessing_log = os.path.join(main_dir, name + "/data/"+name+"_preprocessing.log")
	if os.stat(preprocessing_log).st_size != 0:
			log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".stdout.txt"))

	preprocessing_log = os.path.join(main_dir, name + "/data/"+name+"_DNASE_PE.sh")
	if os.stat(preprocessing_log).st_size != 0:
			log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".script.sh"))

	preprocessing_log = os.path.join(main_dir,  name + "/data/"+name+"_bias_pwm.png")
	if os.stat(preprocessing_log).st_size != 0:
			log_paths.append((preprocessing_log,"logfile.preprocessing."+encid+".bias_pwm.png"))

	# peak-calling-log-files
	tmpdir = "/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/DNASE/caper/"
	
	peaks_log = os.path.join(tmpdir, name + "/metadata.json")
	if os.path.isfile(peaks_log):
			log_paths.append((peaks_log,"logfile.peak_calling."+encid+".metadata.json"))

	peaks_log = os.path.join(tmpdir, name + "/call-reproducibility_overlap/stdout")
	if os.path.isfile(peaks_log):
			log_paths.append((peaks_log,"logfile.peak_calling."+encid+".stdout.txt"))

	peaks_log = os.path.join(tmpdir, name + "/call-reproducibility_overlap/stderr")
	if os.path.isfile(peaks_log):
			log_paths.append((peaks_log,"logfile.peak_calling."+encid+".stderr.txt"))

	return log_paths
	
def bias_fetch_preprocessing_log_files_set_2(odir, encid, main_dir, name):
	# do bed file checks
	log_paths = []
	
	# preprocessing log files
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

	# peak-calling-log-files
	peaks_log = os.path.join(odir, encid + "/peak_calling/log.e")
	if os.path.isfile(peaks_log):
			log_paths.append((peaks_log,"logfile.peak_calling."+encid+".stdout_v1.txt"))

	peaks_log = os.path.join(odir, encid + "/peak_calling/log.o")
	if os.path.isfile(peaks_log):
			log_paths.append((peaks_log,"logfile.peak_calling."+encid+".stdout_v2.txt"))

	return log_paths
		
def fetch_per_fold_training_data_bias(odir, model_dir, encid, fold_num, main_dir, name):
	input_paths = []
	log_paths = []
	
	#print(model_dir)
	opath = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/splits_format/"
	filtered_regions_bed = os.path.join(opath + "/fold_"+str(fold_num)+".json")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"cv_params.fold_"+str(fold_num)+".json"))

	if fold_num==0:
		if name in ["HEPG2"]:
			filtered_regions_bed = os.path.join(main_dir, name+"/negatives_data/negatives_with_summit.bed.gz")
			if os.path.isfile(filtered_regions_bed):
				input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))
		else:

			filtered_regions_bed = os.path.join(main_dir, name+"/data/negatives_data/negatives_with_summit.bed.gz")
			if os.path.isfile(filtered_regions_bed):
				input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))
		
			filtered_regions_bed = os.path.join(odir, encid+"/negatives_data/negatives_with_summit.bed.gz")
			if os.path.isfile(filtered_regions_bed):
				input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	else:
		if name in ["HEPG2"]:
			filtered_regions_bed = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/negatives_with_summit.bed.gz")
			print(filtered_regions_bed)
			if os.path.isfile(filtered_regions_bed):
				input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))

		else:
			filtered_regions_bed = os.path.join(main_dir, name+"/data/negatives_data_"+str(fold_num)+"/negatives_with_summit.bed.gz")
			print(filtered_regions_bed)
			if os.path.isfile(filtered_regions_bed):
				input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))
		

			filtered_regions_bed = os.path.join(odir, encid+"/negatives_data_"+str(fold_num)+"/negatives_with_summit.bed.gz")
			if os.path.isfile(filtered_regions_bed):
				input_paths.append((filtered_regions_bed,"nonpeaks.all_input_regions.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	#filtered_regions_bed = os.path.join(model_dir, "bias_model/train_test_regions/peaks.testset.bed.gz")
	#print(filtered_regions_bed)
	#if os.path.isfile(filtered_regions_bed):
	#	input_paths.append((filtered_regions_bed,"peaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_bias/nonpeaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_bias/nonpeaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions_bias/nonpeaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))
			
	#print(input_paths)
	#print(filtered_regions_bed)
	
	if fold_num==0:
		if name in ["HEPG2"]:
			negatives_log = os.path.join(main_dir, name+"/negatives_data/make_background_regions.log")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))
		else:
	
			negatives_log = os.path.join(main_dir, name+"/data/negatives_data/make_background_regions.log")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))

			negatives_log = os.path.join(odir, encid+"/negatives_data/make_background_regions.log")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))

			negatives_log = os.path.join(odir, encid+"/negatives_data/gc_matching.log.o")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout_v2.txt"))

	else:
		if name in ["HEPG2"]:
			negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/make_background_regions.log")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))
		else:
			negatives_log = os.path.join(main_dir, name+"/data/negatives_data_"+str(fold_num)+"/make_background_regions.log")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))

			negatives_log = os.path.join(odir, encid+"/negatives_data_"+str(fold_num)+"/make_background_regions.log")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout_v1.txt"))


	if fold_num==0:
		if name in ["HEPG2"]:
			negatives_log = os.path.join(main_dir, name+"/negatives_data/negatives_compared_with_foreground.png")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))
		else:
			negatives_log = os.path.join(main_dir, name+"/data/negatives_data/negatives_compared_with_foreground.png")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))

			negatives_log = os.path.join(odir, encid+"/negatives_data/negatives_compared_with_foreground.png")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))

	else:
		if name in ["HEPG2"]:
			negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/negatives_compared_with_foreground.png")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))
		else:
			negatives_log = os.path.join(main_dir, name+"/data/negatives_data_"+str(fold_num)+"/negatives_compared_with_foreground.png")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))

			negatives_log = os.path.join(odir, encid+"/negatives_data_"+str(fold_num)+"/negatives_compared_with_foreground.png")
			if os.path.isfile(negatives_log):
				if os.stat(negatives_log).st_size != 0:
					log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))

#	negatives_log = os.path.join(odir, encid + "/negatives_data/test/fold_"+str(fold_num)+"."+encid+"_test.log")
#	if os.stat(negatives_log).st_size != 0:
#		log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	# add preprocessing data main_dir
	
	return input_paths, log_paths


def fetch_per_fold_bias_models(odir, model_dir, encid, fold_num):
	input_paths = []
	log_paths = []

	bm_model = os.path.join(model_dir, "bias_model/bias.h5")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		return None, None
			
	bm_model = os.path.join(model_dir, "bias_model/new_model_formats_v1/bias.tar")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None
					
	#### fetch model training log files ########
	
	modelling_log = os.path.join(model_dir, "bias_model/train_bias_model.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	modelling_log = os.path.join(model_dir, "bias_model/bias.args.json")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".args.json"))
	
	modelling_log = os.path.join(model_dir, "bias_model/bias_data_params.tsv")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".bias_data_params.tsv"))

	modelling_log = os.path.join(model_dir, "bias_model/bias_model_params.tsv")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".bias_train_params.tsv"))

	modelling_log = os.path.join(model_dir, "bias_model/bias.params.json")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".bias_train_params.json"))

	modelling_log = os.path.join(model_dir, "bias_model/bias.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".epoch_loss.csv"))

	modelling_log = os.path.join(model_dir, "bias_model/bias.log.batch")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile.modelling.fold_"+str(fold_num)+"."+encid+".batch_loss.tsv"))

	return input_paths, log_paths


			
			
	
