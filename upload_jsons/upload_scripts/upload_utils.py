import os
import json
import numpy as np

### utils for preprocessing

		
### utils for training and testing regions

def bias_fetch_preprocessing_log_files(odir, encid, main_dir, name):
	# do bed file checks
	log_paths = []
	# preprocessing, negatives, peak-calling
										
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
	
def fetch_preprocessing_log_files(odir, encid, main_dir, name):
	# do bed file checks
	log_paths = []
	# preprocessing, negatives, peak-calling
	
	temp_dir="/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/"
	preprocessing_log = os.path.join(temp_dir, name + "/script.sh")
	if os.stat(preprocessing_log).st_size != 0:
		log_paths.append((preprocessing_log,"logfile_v1.preprocessing."+encid+".stdout.txt"))

	preprocessing_log = os.path.join(main_dir, name + "/data/"+name+"_preprocessing.log")
	if os.stat(preprocessing_log).st_size != 0:
		log_paths.append((preprocessing_log,"logfile_v2.preprocessing."+encid+".stdout.txt"))
												
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
	
def fetch_per_fold_training_data_bias(odir, model_dir, encid, fold_num, main_dir, name):
	input_paths = []
	log_paths = []
	
	#print(model_dir)
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
		

	filtered_regions_bed = os.path.join(model_dir, "bias_model/train_test_regions/peaks.testset.bed.gz")
	#print(filtered_regions_bed)
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "bias_model/train_test_regions/nonpeaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "bias_model/train_test_regions/nonpeaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "bias_model/train_test_regions/nonpeaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))
			
	#print(input_paths)
	#print(filtered_regions_bed)
	
	if fold_num==0:
		negatives_log = os.path.join(main_dir, name+"/negatives_data/make_background_regions.log")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))
	else:
		negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/make_background_regions.log")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))


	if fold_num==0:
		negatives_log = os.path.join(main_dir, name+"/negatives_data/negatives_compared_with_foreground.png")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))
	else:
		negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/negatives_compared_with_foreground.png")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))

	negatives_log = os.path.join(odir, encid + "/negatives_data/test/fold_"+str(fold_num)+"."+encid+"_test.log")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	# add preprocessing data main_dir
	
	return input_paths, log_paths

def fetch_per_fold_training_data(odir, model_dir, encid, fold_num, main_dir, name):
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
		
	filtered_regions_bed = os.path.join(model_dir, "train_test_regions/peaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))


	filtered_regions_bed = os.path.join(model_dir, "train_test_regions/peaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions/peaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"peaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions/nonpeaks.trainingset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.trainingset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions/nonpeaks.validationset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.validationset.fold_"+str(fold_num)+"."+encid+".bed.gz"))

	filtered_regions_bed = os.path.join(model_dir, "train_test_regions/nonpeaks.testset.bed.gz")
	if os.path.isfile(filtered_regions_bed):
		input_paths.append((filtered_regions_bed,"nonpeaks.testset.fold_"+str(fold_num)+"."+encid+".bed.gz"))
			
	#input_paths)
	#print(filtered_regions_bed)
	
	if fold_num==0:
		negatives_log = os.path.join(main_dir, name+"/negatives_data/make_background_regions.log")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))
	else:
		negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/make_background_regions.log")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))


	if fold_num==0:
		negatives_log = os.path.join(main_dir, name+"/negatives_data/negatives_compared_with_foreground.png")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))
	else:
		negatives_log = os.path.join(main_dir, name+"/negatives_data_"+str(fold_num)+"/negatives_compared_with_foreground.png")
		if os.stat(negatives_log).st_size != 0:
			log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.png"))

	negatives_log = os.path.join(odir, encid + "/negatives_data/test/fold_"+str(fold_num)+"."+encid+"_test.log")
	if os.stat(negatives_log).st_size != 0:
		log_paths.append((negatives_log,"logfile.gc_matching.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	# add preprocessing data main_dir
	
	return input_paths, log_paths

### utils for model uploads

#just need to add log files

def fetch_per_fold_bias_models(odir, model_dir, encid, fold_num):
	input_paths = []
	log_paths = []

	bm_model = os.path.join(model_dir, "bias_model/bias.h5")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		return None, None
			
	bm_model = os.path.join(model_dir, "bias_model/new_model_formats/bias.tar")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None
					
	#### fetch model training log files ########

	modelling_log = os.path.join(model_dir, "bias_model/bias_data_params.tsv")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile_v1.modelling.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	modelling_log = os.path.join(model_dir, "bias_model/bias.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile_v2.modelling.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	return input_paths, log_paths

def fetch_per_fold_models(odir, model_dir, encid, fold_num):
	input_paths = []
	log_paths = []
	
	cmb = os.path.join(model_dir, "chrombpnet_model/chrombpnet_wo_bias.h5")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet_nobias.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		return None, None
		
	cm_model = os.path.join(model_dir, "chrombpnet_model/chrombpnet.h5")
	if os.path.isfile(cm_model):
		input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		return None, None
			
	bm_model = os.path.join(model_dir, "chrombpnet_model/bias_model_scaled.h5")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias_scaled.fold_"+str(fold_num)+"."+encid+".h5"))
	else:
		return None, None

	cmb = os.path.join(model_dir, "chrombpnet_model/new_model_formats/chrombpnet_wo_bias.tar")
	if os.path.isfile(cmb):
		input_paths.append((cmb,"model.chrombpnet_nobias.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None
		
	cm_model = os.path.join(model_dir, "chrombpnet_model/new_model_formats/chrombpnet.tar")
	if os.path.isfile(cm_model):
		input_paths.append((cm_model,"model.chrombpnet.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None
			
	bm_model = os.path.join(model_dir, "chrombpnet_model/new_model_formats/bias_model_scaled.tar")
	if os.path.isfile(bm_model):
		input_paths.append((bm_model,"model.bias_scaled.fold_"+str(fold_num)+"."+encid+".tar"))
	else:
		return None, None
					
	#### fetch model training log files ########

	modelling_log = os.path.join(model_dir, "chrombpnet_model/train_chrombpnet_model.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile_v1.modelling.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	modelling_log = os.path.join(model_dir, "chrombpnet_model/chrombpnet.log")
	if os.stat(modelling_log).st_size != 0:
		log_paths.append((modelling_log,"logfile_v2.modelling.fold_"+str(fold_num)+"."+encid+".stdout.txt"))

	return input_paths, log_paths

### utils for model performance metrics uploads

def fetch_chrombpnet_median_performance(odir, encid, args_json, model_dirs, main_dir, name):

	success = False
	
	os.makedirs(os.path.join(main_dir,name+"/median_performance/"), exist_ok=True)
	
	cmb = os.path.join(main_dir, name + "/median_performance/median_chrombpnet_performance.json")
	if os.path.isfile(cmb):
		success = True
		return success, cmb	 
	else:
		median_metrics = {}
		median_metrics["counts_metrics"] = {"peaks": {}, "nonpeaks": {}, "peaks_and_nonpeaks":{}}
		median_metrics["profile_metrics"] = {"peaks": {}, "nonpeaks": {}, "peaks_and_nonpeaks":{}}
		
		summary = {}
		for fold_num in range(5):
			json_path = os.path.join(model_dirs[fold_num], "chrombpnet_model/new_test_metrics/chrombpnet_metrics.json")
			if os.path.isfile(json_path):
				metrics_dict = json.load(open(json_path))
				for key1 in metrics_dict.keys():
					for key2 in metrics_dict[key1].keys():
						for key3 in metrics_dict[key1][key2].keys():
							if key1+"."+key2+"."+key3 not in summary:
								summary[key1+"."+key2+"."+key3] = []
							summary[key1+"."+key2+"."+key3].append(metrics_dict[key1][key2][key3])
			
		if len(summary) == 0:
			success = False
			return success, cmb
			
		for allk in summary:
			#print(len(summary[allk]))
			if len(summary[allk]) != 5:
				success = False
				return success, cmb
				
			val = np.median(summary[allk])
			keysp = allk.split(".")
			median_metrics[keysp[0]][keysp[1]][keysp[2]] = val
			
					
		with open(cmb, "w") as ofile:
			json.dump(median_metrics, ofile, indent=4)
		
		success = True
		return success, cmb	

def fetch_chrombpnet_nobias_median_performance(odir, encid, args_json, model_dirs, main_dir, name):

	success = False
	
	os.makedirs(os.path.join(main_dir,name+"/median_performance/"), exist_ok=True)
	
	cmb = os.path.join(main_dir, name + "/median_performance/median_chrombpnet_performance.json")
	if os.path.isfile(cmb):
		success = True
		return success, cmb
	else:
		median_metrics = {}
		median_metrics["counts_metrics"] = {"peaks": {}, "nonpeaks": {}, "peaks_and_nonpeaks":{}}
		median_metrics["profile_metrics"] = {"peaks": {}, "nonpeaks": {}, "peaks_and_nonpeaks":{}}
		
		summary = {}
		for fold_num in range(5):
			json_path = os.path.join(odir, encid + "/" + model_dirs[fold_num] + "/new_test_metrics/chrombpnet_wo_bias_metrics.json")
			if os.path.isfile(json_path):
				metrics_dict = json.load(open(json_path))
				for key1 in metrics_dict.keys():
					for key2 in metrics_dict[key1].keys():
						for key3 in metrics_dict[key1][key2].keys():
							if key1+"."+key2+"."+key3 not in summary:
								summary[key1+"."+key2+"."+key3] = []
							summary[key1+"."+key2+"."+key3].append(metrics_dict[key1][key2][key3])

		if len(summary) == 0:
			success = False
			return success, cmb
						
		for allk in summary:
			#print(len(summary[allk]))
			if len(summary[allk]) != 5:
				success = False
				return success, cmb
				
			val = np.median(summary[allk])
			keysp = allk.split(".")
			median_metrics[keysp[0]][keysp[1]][keysp[2]] = val
						
		with open(cmb, "w") as ofile:
			json.dump(median_metrics, ofile, indent=4)
		
		success = True
		return success, cmb

def fetch_model_perf_logs(odir, encid, name):
	log_paths = []

	temp_dir="/oak/stanford/groups/akundaje/projects/chrombpnet/model_inputs/ENCODE_ATAC_downloads/"
	preprocessing_log = os.path.join(temp_dir, name + "/samstats.sh")
	if os.path.isfile(preprocessing_log):
		log_paths.append((preprocessing_log,"logfile.samstats."+encid+".stdout.txt"))

	return log_paths
	

def fetch_per_fold_model_metrics(odir, model_dir, encid, fold_num, ddtype):
	data_paths = []
	log_paths = []

	# predictions h5
	
	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_predictions.h5")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "predictions.chrombpnet.testset.fold_"+str(fold_num)+"."+encid+".h5"))

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_wo_bias_predictions.h5")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "predictions.chrombpnet_nobias.testset.fold_"+str(fold_num)+"."+encid+".h5"))

	# only peaks scatter
	
	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_only_peaks.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "obs_pred_scatter_plot.counts.chrombpnet.peaks.testset.fold_"+str(fold_num)+"."+encid+".png"))

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_wo_bias_only_peaks.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "obs_pred_scatter_plot.counts.chrombpnet_nobias.peaks.testset.fold_"+str(fold_num)+"."+encid+".png"))

	# all regions scatter
	
	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_peaks_and_nonpeaks.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "obs_pred_scatter_plot.counts.chrombpnet.all.testset.fold_"+str(fold_num)+"."+encid+".png"))

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_wo_bias_peaks_and_nonpeaks.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "obs_pred_scatter_plot.counts.chrombpnet_nobias.all.testset.fold_"+str(fold_num)+"."+encid+".png"))

	# only peaks jsd

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_only_peaks.jsd.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "jsd_histogram_plot.profiles.chrombpnet.peaks.testset.fold_"+str(fold_num)+"."+encid+".png"))

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_wo_bias_only_peaks.jsd.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "jsd_histogram_plot.profiles.chrombpnet_nobias.peaks.testset.fold_"+str(fold_num)+"."+encid+".png"))

	# all regions jsd

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_peaks_and_nonpeaks.jsd.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "jsd_histogram_plot.profiles.chrombpnet.all.testset.fold_"+str(fold_num)+"."+encid+".png"))

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_wo_bias_peaks_and_nonpeaks.jsd.png")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "jsd_histogram_plot.profiles.chrombpnet_nobias.all.testset.fold_"+str(fold_num)+"."+encid+".png"))
	
	# summary json

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_metrics.json")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "performance_metrics.chrombpnet.fold_"+str(fold_num)+"."+encid+".json"))

	cmb = os.path.join(model_dir, "chrombpnet_model/new_test_metrics/chrombpnet_wo_bias_metrics.json")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "performance_metrics.chrombpnet_nobias.fold_"+str(fold_num)+"."+encid+".json"))

	
	# footprint summary
	
	cmb = os.path.join(model_dir, "chrombpnet_model/footprints/corrected_footprints_score.txt")
	if os.path.isfile(cmb):
		data_paths.append((cmb, "footprint_plot.chrombpnet_nobias.bias_motif.summary.fold_"+str(fold_num)+"."+encid+".txt"))
	else:
		cmb = os.path.join(model_dir, "chrombpnet_model/footprints_motifs/motif_footprints_score.txt")
		if os.path.isfile(cmb):
			data_paths.append((cmb, "footprint_plot.chrombpnet_nobias.bias_motif.summary.fold_"+str(fold_num)+"."+encid+".txt"))

	# footprint folder files
	
	footprint_paths = []
	
	if ddtype=="ATAC":
		for j in range(1,6):
			cmb = os.path.join(model_dir, "chrombpnet_model/footprints/corrected.tn5_"+str(j)+".footprint.png")
			if os.path.isfile(cmb):
				footprint_paths.append((cmb, "footprint.chrombpnet_nobias.bias_motif_"+str(j)+".fold_"+str(fold_num)+"."+encid+".png"))
			else:
				cmb = os.path.join(model_dir, "chrombpnet_model/footprints_motifs/motif.tn5_"+str(j)+".footprint.png")
				if os.path.isfile(cmb):
					footprint_paths.append((cmb, "footprint.chrombpnet_nobias.bias_motif_"+str(j)+".fold_"+str(fold_num)+"."+encid+".png"))
		
	
			#cmb = os.path.join(model_dir, "bias_model/footprints/bias.tn5_"+str(j)+".footprint.png")
			#if os.path.isfile(cmb):
			#	footprint_paths.append((cmb, "footprint.bias_scaled.bias_motif_"+str(j)+".fold_"+str(fold_num)+"."+encid+".png"))

	if ddtype=="DNASE":
		for j in range(1,3):
			cmb = os.path.join(model_dir, "chrombpnet_model/footprints/corrected.dnase_"+str(j)+".footprint.png")
			if os.path.isfile(cmb):
				footprint_paths.append((cmb, "footprint.chrombpnet_nobias.bias_motif_"+str(j)+".fold_"+str(fold_num)+"."+encid+".png"))
			else:
				cmb = os.path.join(model_dir, "chrombpnet_model/footprints_motifs/motif.dnase_"+str(j)+".footprint.png")
				if os.path.isfile(cmb):
					footprint_paths.append((cmb, "footprint.chrombpnet_nobias.bias_motif_"+str(j)+".fold_"+str(fold_num)+"."+encid+".png"))
		
			#cmb = os.path.join(model_dir, "/footprints/bias.dnase_"+str(j)+".footprint.png")
			#if os.path.isfile(cmb):
			#	footprint_paths.append((cmb, "footprint.bias_scaled.bias_motif_"+str(j)+".fold_"+str(fold_num)+"."+encid+".png"))
	
	#cmb = os.path.join(model_dir, "/footprints/bias_footprints_score.txt")
	#if os.path.isfile(cmb):
	#	footprint_paths.append((cmb, "footprint_plot.bias_scaled.bias_motif.summary.fold_"+str(fold_num)+"."+encid+".txt"))

	#cmb = os.path.join(model_dir, "/footprints/bias_footprints.h5")
	#if os.path.isfile(cmb):
#		footprint_paths.append((cmb, "footprint.bias_scaled.bias_motif.all.fold_"+str(fold_num)+"."+encid+".h5"))
	
	cmb = os.path.join(model_dir, "chrombpnet_model/footprints/corrected_footprints.h5")
	if os.path.isfile(cmb):
		footprint_paths.append((cmb, "footprint.chrombpnet_nobias.bias_motif.all.fold_"+str(fold_num)+"."+encid+".h5"))
	else:		
		cmb = os.path.join(model_dir, "chrombpnet_model/footprints_motifs/motif_footprints.h5")
		if os.path.isfile(cmb):
			footprint_paths.append((cmb, "footprint.chrombpnet_nobias.bias_motif.all.fold_"+str(fold_num)+"."+encid+".h5"))
		
	return 	data_paths, log_paths, footprint_paths

def fetch_mean_metrics_json(odir, models_path, encid):
	mean_metrics = {}
	for model_dir in models_path:
		cmb = os.path.join(odir, encid + "/" + model_dir + "/test/chrombpnet_metrics.json")
		metrics = json.load(open(cmb))
		for key in metrics:
			if key not in mean_metrics:
				mean_metrics[key] = metrics[key]
			else:
				mean_metrics[key] += metrics[key]	
			
	for key in mean_metrics:
		mean_metrics[key] = mean_metrics[key] / 5
	
	output_path = os.path.join(odir, encid + "/" + "bpnet_quality_metrics.json")
	with open(output_path, "w") as outfile:
		json.dump(mean_metrics, outfile, indent=4)

	return output_path


			
			
	
