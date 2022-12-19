import h5py
import hdf5plugin
import argparse
import chrombpnet.parsers as parsers
import numpy as np
import pandas as pd
import os
import json
import copy
from chrombpnet.data import DefaultDataFile, get_default_data_path
from chrombpnet.data import print_meme_motif_file

def chrombpnet_train_pipeline(args):

	# Shift bam and convert to bigwig
	import chrombpnet.helpers.preprocessing.reads_to_bigwig as reads_to_bigwig	
	args.output_prefix = os.path.join(args.output_dir,"auxiliary/data")
	reads_to_bigwig.main(args)
	
	# QC bigwig
	import chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig as build_pwm_from_bigwig	
	args.bigwig = os.path.join(args.output_dir,"auxiliary/data_unstranded.bw")
	args.output_prefix = os.path.join(args.output_dir,"evaluation/bw_shift_qc")
	folds = json.load(open(args.chr_fold_path))
	assert(len(folds["valid"]) > 0) # validation list of chromosomes is empty
	args.chr = folds["valid"][0]
	args.pwm_width=24
	build_pwm_from_bigwig.main(args)
	
	# make predictions with input bias model in peaks
	import chrombpnet.training.predict as predict
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args_copy.output_dir,"evaluation/bias")
	args_copy.model_h5 = args.bias_model_path
	args_copy.nonpeaks = "None"
	predict.main(args_copy)
	
	# QC bias model performance in peaks
	bias_metrics = json.load(open(os.path.join(args_copy.output_dir,"evaluation/bias_metrics.json")))
	assert(bias_metrics["counts_metrics"]["peaks"]["pearsonr"] > -0.3) # bias model has negative correlation in peaks - AT rich bias model. Increase bias threshold and retrain bias model. Or use a different bias model with higher bias threshold. 
	
	# fetch hyperparameters for training
	import chrombpnet.helpers.hyperparameters.find_chrombpnet_hyperparams as find_chrombpnet_hyperparams
	args_copy = copy.deepcopy(args)
	args_copy.output_dir = os.path.join(args.output_dir,"auxiliary/")
	find_chrombpnet_hyperparams.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"auxiliary/bias_model_scaled.h5"),os.path.join(args.output_dir,"models/bias_model_scaled.h5"))
	os.rename(os.path.join(args.output_dir,"auxiliary/chrombpnet_model_params.tsv"),os.path.join(args.output_dir,"logs/chrombpnet_model_params.tsv"))
	os.rename(os.path.join(args.output_dir,"auxiliary/chrombpnet_data_params.tsv"),os.path.join(args.output_dir,"logs/chrombpnet_data_params.tsv"))

	params = open(os.path.join(args.output_dir,"logs/chrombpnet_model_params.tsv")).read()
	params = params.replace(os.path.join(args.output_dir,"auxiliary/bias_model_scaled.h5"),os.path.join(args.output_dir,"models/bias_model_scaled.h5"))
	with open(os.path.join(args.output_dir,"logs/chrombpnet_model_params.tsv"),"w") as f:
		f.write(params)
		
	# get model architecture path and train chromBPNet model
	import chrombpnet.training.models.chrombpnet_with_bias_model as chrombpnet_with_bias_model
	import chrombpnet.training.train as train
	args_copy = copy.deepcopy(args)
	if args_copy.architecture_from_file is None:
		args_copy.architecture_from_file = 	chrombpnet_with_bias_model.__file__
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/filtered.peaks.bed")
	args_copy.nonpeaks = os.path.join(args.output_dir,"auxiliary/filtered.nonpeaks.bed")
	args_copy.output_prefix = os.path.join(args.output_dir,"models/chrombpnet")
	args_copy.params = os.path.join(args.output_dir,"logs/chrombpnet_model_params.tsv")
	train.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"models/chrombpnet.log"),os.path.join(args.output_dir,"logs/chrombpnet.log"))
	os.rename(os.path.join(args.output_dir,"models/chrombpnet.log.batch"),os.path.join(args.output_dir,"logs/chrombpnet.log.batch"))
	os.rename(os.path.join(args.output_dir,"models/chrombpnet.params.json"),os.path.join(args.output_dir,"logs/chrombpnet.params.json"))
	os.rename(os.path.join(args.output_dir,"models/chrombpnet.args.json"),os.path.join(args.output_dir,"logs/chrombpnet.args.json"))

	if args.cmd == "train":
		print("Finished training! Exiting!")
		return
		
	# make predictions with trained chrombpnet model
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args.output_dir,"evaluation/chrombpnet")
	args_copy.model_h5 = os.path.join(args.output_dir,"models/chrombpnet_nobias.h5")
	args_copy.nonpeaks = "None"
	predict.main(args_copy)
	
	# marginal footprinting with model
	import chrombpnet.evaluation.marginal_footprints.marginal_footprinting as marginal_footprinting
	if args.data_type == "ATAC":
		bias_motifs = [["tn5_1","GCACAGTACAGAGCTG"],["tn5_2","GTGCACAGTTCTAGAGTGTGCAG"],["tn5_3","CCTCTACACTGTGCAGAA"],["tn5_4","GCACAGTTCTAGACTGTGCAG"],["tn5_5","CTGCACAGTGTAGAGTTGTGC"]]

	elif args.data_type == "DNASE":
		bias_motifs = [["dnase_1","TTTACAAGTCCA"],["dnase_2","TTTACAAGTCCA"]]
	else:
		print("unknown data type: "+args.data_type)
	df = pd.DataFrame(bias_motifs)
	df.to_csv(os.path.join(args_copy.output_dir,"auxiliary/motif_to_pwm.tsv"),sep="\t",header=False,index=False)
	
	args_copy = copy.deepcopy(args)
	args_copy.model_h5 = os.path.join(args.output_dir,"models/chrombpnet_nobias.h5")
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/filtered.nonpeaks.bed")
	args_copy.output_prefix = os.path.join(args.output_dir,"evaluation/chrombpnet_nobias")
	args_copy.motifs_to_pwm = os.path.join(args_copy.output_dir,"auxiliary/motif_to_pwm.tsv")
	args_copy.ylim = None
	marginal_footprinting.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"evaluation/chrombpnet_nobias_footprints.h5"),os.path.join(args.output_dir,"auxiliary/chrombpnet_nobias_footprints.h5"))

	# get contributions scores with model
	args_copy = copy.deepcopy(args)
	import chrombpnet.evaluation.interpret.interpret as interpret
	peaks = pd.read_csv(os.path.join(args_copy.output_dir,"auxiliary/filtered.peaks.bed"),sep="\t",header=None)
	if peaks.shape[0] > 30000:
		sub_peaks = peaks.sample(30000, random_state=1234)
	else:
		sub_peaks = peaks.sample(30000, random_state=1234)
	sub_peaks.to_csv(os.path.join(args_copy.output_dir,"auxiliary/30K_subsample_peaks.bed"),sep="\t", header=False, index=False)
	os.makedirs(os.path.join(args.output_dir,"auxiliary/interpret/"), exist_ok=False)

	args_copy.profile_or_counts = ["counts", "profile"]
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/30K_subsample_peaks.bed")	
	args_copy.model_h5 = os.path.join(args.output_dir,"models/chrombpnet_nobias.h5")
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/interpret/chrombpnet_nobias")
	args_copy.debug_chr = None
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
	# modisco-lite pipeline
	
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret/chrombpnet_nobias.profile_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_profile_scores.h5"))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_profile_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_profile/"),meme_file)
	os.system(modisco_command)
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret/chrombpnet_nobias.counts_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_counts_scores.h5"))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_counts_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_counts/"),meme_file)
	os.system(modisco_command)

	import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
	convert_html_to_pdf.main(os.path.join(args.output_dir,"auxiliary/interpret/modisco_counts/motifs.html"),os.path.join(args.output_dir,"evaluation/chrombpnet_nobias_counts.pdf"))
	convert_html_to_pdf.main(os.path.join(args.output_dir,"auxiliary/interpret/modisco_profile/motifs.html"),os.path.join(args.output_dir,"evaluation/chrombpnet_nobias_profile.pdf"))
	
	import chrombpnet.helpers.generate_reports.make_html as make_html
	args_copy = copy.deepcopy(args)
	args_copy.input_dir = args_copy.output_dir
	make_html(args_copy)

def train_bias_pipeline(args):

	# Shift bam and convert to bigwig
	import chrombpnet.helpers.preprocessing.reads_to_bigwig as reads_to_bigwig	
	args.output_prefix = os.path.join(args.output_dir,"auxiliary/data")
	reads_to_bigwig.main(args)
	
	# QC bigwig
	import chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig as build_pwm_from_bigwig	
	args.bigwig = os.path.join(args.output_dir,"auxiliary/data_unstranded.bw")
	args.output_prefix = os.path.join(args.output_dir,"evaluation/bw_shift_qc")
	folds = json.load(open(args.chr_fold_path))
	assert(len(folds["valid"]) > 0) # validation list of chromosomes is empty
	args.chr = folds["valid"][0]
	args.pwm_width=24
	build_pwm_from_bigwig.main(args)
	

	# fetch hyperparameters for training
	import chrombpnet.helpers.hyperparameters.find_bias_hyperparams as find_bias_hyperparams
	args_copy = copy.deepcopy(args)
	args_copy.output_dir = os.path.join(args.output_dir,"auxiliary/")
	find_bias_hyperparams.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"auxiliary/bias_model_params.tsv"),os.path.join(args.output_dir,"logs/bias_model_params.tsv"))
	os.rename(os.path.join(args.output_dir,"auxiliary/bias_data_params.tsv"),os.path.join(args.output_dir,"logs/bias_data_params.tsv"))
		
	# get model architecture path and train chromBPNet model
	import chrombpnet.training.models.bpnet_model as bpnet_model
	import chrombpnet.training.train as train
	args_copy = copy.deepcopy(args)
	if args_copy.architecture_from_file is None:
		args_copy.architecture_from_file = 	bpnet_model.__file__
	args_copy.peaks = "None"
	args_copy.nonpeaks = os.path.join(args_copy.output_dir,"auxiliary/filtered.bias_nonpeaks.bed")
	args_copy.output_prefix = os.path.join(args_copy.output_dir,"models/bias")
	args_copy.params = os.path.join(args_copy.output_dir,"logs/bias_model_params.tsv")
	train.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"models/bias.args.json"),os.path.join(args.output_dir,"logs/bias.args.json"))
	os.rename(os.path.join(args.output_dir,"models/bias.log"),os.path.join(args.output_dir,"logs/bias.log"))
	os.rename(os.path.join(args.output_dir,"models/bias.log.batch"),os.path.join(args.output_dir,"logs/bias.log.batch"))
	os.rename(os.path.join(args.output_dir,"models/bias.params.json"),os.path.join(args.output_dir,"logs/bias.params.json"))

	if args.cmd_bias == "train":
		print("Finished training! Exiting!")
		return
		
	# make predictions with trained bias model
	import chrombpnet.training.predict as predict
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args_copy.output_dir,"evaluation/bias")
	args_copy.model_h5 = os.path.join(args.output_dir,"models/bias.h5")
	args_copy.peaks = "None"
	predict.main(args_copy)
	
	# get contributions scores with model
	import chrombpnet.evaluation.interpret.interpret as interpret
	peaks = pd.read_csv(os.path.join(args.peaks),sep="\t",header=None)
	if peaks.shape[0] > 30000:
		sub_peaks = peaks.sample(30000, random_state=1234)
	else:
		sub_peaks = peaks.sample(30000, random_state=1234)
	sub_peaks.to_csv(os.path.join(args_copy.output_dir,"auxiliary/30K_subsample_peaks.bed"),sep="\t", header=False, index=False)
	os.makedirs(os.path.join(args.output_dir,"auxiliary/interpret/"), exist_ok=False)

	args_copy = copy.deepcopy(args)
	args_copy.profile_or_counts = ["counts", "profile"]
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/30K_subsample_peaks.bed")	
	args_copy.model_h5 = os.path.join(args.output_dir,"models/bias.h5")
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/interpret/bias")
	args_copy.debug_chr = None
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	# modisco-lite pipeline
	
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret/bias.profile_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_profile_scores.h5"))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_profile_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_profile/"),meme_file)
	os.system(modisco_command)
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret/bias.counts_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_counts_scores.h5"))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret/modisco_results_counts_scores.h5"),os.path.join(args.output_dir,"auxiliary/interpret/modisco_counts/"),meme_file)
	os.system(modisco_command)
	
	import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
	convert_html_to_pdf.main(os.path.join(args.output_dir,"auxiliary/interpret/modisco_counts/motifs.html"),os.path.join(args.output_dir,"evaluation/bias_counts.pdf"))
	convert_html_to_pdf.main(os.path.join(args.output_dir,"auxiliary/interpret/modisco_profile/motifs.html"),os.path.join(args.output_dir,"evaluation/bias_profile.pdf"))

	import chrombpnet.helpers.generate_reports.make_html_bias as make_html_bias
	args_copy = copy.deepcopy(args)
	args_copy.input_dir = args_copy.output_dir
	make_html(args_copy_bias)	
		
def main():
	args = parsers.read_parser()
	if args.cmd == "pipeline" or args.cmd == "train":
		os.makedirs(os.path.join(args.output_dir,"logs"), exist_ok=False)
		os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=False)
		os.makedirs(os.path.join(args.output_dir,"models"), exist_ok=False)
		os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=False)

		chrombpnet_train_pipeline(args)
		
	elif args.cmd == "bias":
		if args.cmd_bias == "pipeline" or args.cmd_bias == "train":
			os.makedirs(os.path.join(args.output_dir,"logs"), exist_ok=False)
			os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=False)
			os.makedirs(os.path.join(args.output_dir,"models"), exist_ok=False)
			os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=False)

			train_bias_pipeline(args)
	
	elif args.cmd == "pred_bw":
	
		assert (args.bias_model is None) + (args.chrombpnet_model is None) + (args.chrombpnet_model_nb is None) < 3, "No input model provided!"
		import chrombpnet.evaluation.make_bigwigs.predict_to_bigwig as predict_to_bigwig

		predict_to_bigwig.main(args)

	elif args.cmd == "contribs_bw":
	
		import chrombpnet.evaluation.interpret.interpret as interpret
		interpret.main(args)
		import chrombpnet.evaluation.make_bigwigs.importance_hdf5_to_bigwig as importance_hdf5_to_bigwig
		if "counts" in  args.profile_or_counts:
			args_copy = copy.deepcopy(args)
			args_copy.hdf5 = args_copy.output_prefix + ".counts_scores.h5"
			args_copy.output_prefix = args.output_prefix + ".counts_scores"
			
			importance_hdf5_to_bigwig.main(args_copy)
		if "profile" in  args.profile_or_counts:
			args_copy = copy.deepcopy(args)
			args_copy.hdf5 = args_copy.output_prefix + ".profile_scores.h5"
			args_copy.output_prefix = args.output_prefix + ".profile_scores"
	
			importance_hdf5_to_bigwig.main(args_copy)
			
	elif args.cmd == "footprints":
	
		import chrombpnet.evaluation.marginal_footprints.marginal_footprinting as marginal_footprinting
		marginal_footprinting.main(args)

	elif args.cmd == "snp_score":
	
		import chrombpnet.evaluation.variant_effect_prediction.snp_scoring as snp_scoring
		snp_scoring.main(args)
		
	elif args.cmd == "modisco_motifs":
		import chrombpnet
		chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
		meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
		modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(args.h5py,args.output_prefix+"_modisco.h5")
		os.system(modisco_command)
		modisco_command = "modisco report -i {} -o {} -m {}".format(args.output_prefix+"_modisco.h5",args.output_prefix+"_reports",meme_dir)
		os.system(modisco_command)
		
		import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
		convert_html_to_pdf.main(args.output_prefix+"_reports/motifs.html",args.output_prefix+"_reports/motifs.pdf")

		
	else:
		print("Command not found")


if __name__=="__main_-":
	main()

    
        
