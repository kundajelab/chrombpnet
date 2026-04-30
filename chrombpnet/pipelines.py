import pandas as pd
import os
import json
import copy
from chrombpnet.data import DefaultDataFile, get_default_data_path
from chrombpnet.data import print_meme_motif_file
import numpy as np

def chrombpnet_train_pipeline(args):

	if args.file_prefix:
		fpx = args.file_prefix+"_"
	else:
		fpx = ""
		
	print(args)
	# Shift bam and convert to bigwig
	import chrombpnet.helpers.preprocessing.reads_to_bigwig as reads_to_bigwig	
	args.output_prefix = os.path.join(args.output_dir,"auxiliary/{}data".format(fpx))
	# args.plus_shift = None
	# args.minus_shift = None
	reads_to_bigwig.main(args)
	
	# QC bigwig
	import chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig as build_pwm_from_bigwig	
	args.bigwig = os.path.join(args.output_dir,"auxiliary/{}data_unstranded.bw".format(fpx))
	args.output_prefix = os.path.join(args.output_dir,"evaluation/{}bw_shift_qc".format(fpx))
	folds = json.load(open(args.chr_fold_path))
	assert(len(folds["valid"]) > 0) # validation list of chromosomes is empty
	args.chr = folds["valid"][0]
	args.pwm_width=24
	build_pwm_from_bigwig.main(args)

	# fetch hyperparameters for training
	import chrombpnet.helpers.hyperparameters.find_chrombpnet_hyperparams as find_chrombpnet_hyperparams
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/{}".format(fpx))
	find_chrombpnet_hyperparams.main(args_copy)
	
	# make predictions with input bias model in peaks
	import chrombpnet.training.predict as predict
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args_copy.output_dir,"evaluation/bias")
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/{}filtered.peaks.bed".format(fpx))
	args_copy.model_h5 = args.bias_model_path
	args_copy.nonpeaks = "None"
	predict.main(args_copy)
	
	# QC bias model performance in peaks
	bias_metrics = json.load(open(os.path.join(args_copy.output_dir,"evaluation/bias_metrics.json")))
	print("Bias model pearsonr performance in peaks is: {}".format(str(np.round(bias_metrics["counts_metrics"]["peaks"]["pearsonr"],2))))
	assert(bias_metrics["counts_metrics"]["peaks"]["pearsonr"] > -0.5) # bias model has negative correlation in peaks - AT rich bias model. Increase bias threshold and retrain bias model. Or use a different bias model with higher bias threshold. 
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"auxiliary/{}bias_model_scaled.h5".format(fpx)),os.path.join(args.output_dir,"models/{}bias_model_scaled.h5".format(fpx)))
	os.rename(os.path.join(args.output_dir,"auxiliary/{}chrombpnet_model_params.tsv".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet_model_params.tsv".format(fpx)))
	os.rename(os.path.join(args.output_dir,"auxiliary/{}chrombpnet_data_params.tsv".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet_data_params.tsv".format(fpx)))

	params = open(os.path.join(args.output_dir,"logs/{}chrombpnet_model_params.tsv".format(fpx))).read()
	params = params.replace(os.path.join(args.output_dir,"auxiliary/{}bias_model_scaled.h5".format(fpx)),os.path.join(args.output_dir,"models/{}bias_model_scaled.h5".format(fpx)))
	with open(os.path.join(args.output_dir,"logs/{}chrombpnet_model_params.tsv".format(fpx)),"w") as f:
		f.write(params)
		
	# get model architecture path and train chromBPNet model
	import chrombpnet.training.models.chrombpnet_with_bias_model as chrombpnet_with_bias_model
	import chrombpnet.training.train as train
	args_copy = copy.deepcopy(args)
	if args_copy.architecture_from_file is None:
		args_copy.architecture_from_file = 	chrombpnet_with_bias_model.__file__
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/{}filtered.peaks.bed".format(fpx))
	args_copy.nonpeaks = os.path.join(args.output_dir,"auxiliary/{}filtered.nonpeaks.bed".format(fpx))
	args_copy.output_prefix = os.path.join(args.output_dir,"models/{}chrombpnet".format(fpx))
	args_copy.params = os.path.join(args.output_dir,"logs/{}chrombpnet_model_params.tsv".format(fpx))
	print("args_copy:",args_copy)
	train.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"models/{}chrombpnet.log".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet.log".format(fpx)))
	os.rename(os.path.join(args.output_dir,"models/{}chrombpnet.log.batch".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet.log.batch".format(fpx)))
	os.rename(os.path.join(args.output_dir,"models/{}chrombpnet.params.json".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet.params.json".format(fpx)))
	os.rename(os.path.join(args.output_dir,"models/{}chrombpnet.args.json".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet.args.json".format(fpx)))

	if args.cmd == "train":
		import chrombpnet.helpers.generate_reports.make_html as make_html
		args_copy = copy.deepcopy(args)
		args_copy.input_dir = args_copy.output_dir
		args_copy.command = args.cmd
		make_html.main(args_copy)
		print("Finished training! Exiting!")
		return
		
	# make predictions with trained chrombpnet model
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args.output_dir,"evaluation/{}chrombpnet".format(fpx))
	args_copy.model_h5 = os.path.join(args.output_dir,"models/{}chrombpnet.h5".format(fpx))
	args_copy.nonpeaks = "None"
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/{}filtered.peaks.bed".format(fpx))
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
	args_copy.model_h5 = os.path.join(args.output_dir,"models/{}chrombpnet_nobias.h5".format(fpx))
	args_copy.regions = args.nonpeaks
	args_copy.output_prefix = os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias".format(fpx))
	args_copy.motifs_to_pwm = os.path.join(args_copy.output_dir,"auxiliary/motif_to_pwm.tsv")
	args_copy.ylim = None
	marginal_footprinting.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias_footprints.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/{}chrombpnet_nobias_footprints.h5".format(fpx)))

	# get contributions scores with model
	args_copy = copy.deepcopy(args)
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/{}filtered.peaks.bed".format(fpx))
	import chrombpnet.evaluation.interpret.interpret as interpret
	peaks = pd.read_csv(os.path.join(args_copy.peaks),sep="\t",header=None)
	if peaks.shape[0] > 30000:
		sub_peaks = peaks.sample(30000, random_state=1234)
	else:
		sub_peaks = peaks
	sub_peaks.to_csv(os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx)),sep="\t", header=False, index=False)
	os.makedirs(os.path.join(args.output_dir,"auxiliary/interpret_subsample/"), exist_ok=False)

	#args_copy.profile_or_counts = ["counts", "profile"]
	args_copy.profile_or_counts = ["counts", "profile"]	
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx))	
	args_copy.model_h5 = os.path.join(args.output_dir,"models/{}chrombpnet_nobias.h5".format(fpx))
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias".format(fpx))
	args_copy.debug_chr = None
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
	# modisco-lite pipeline
	
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias.profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_profile/"),meme_file)
	os.system(modisco_command)
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias.counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_counts/"),meme_file)
	os.system(modisco_command)

	import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_counts/motifs.html"),os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias_counts.pdf".format(fpx)))
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_profile/motifs.html"),os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias_profile.pdf".format(fpx)))
	
	import chrombpnet.helpers.generate_reports.make_html as make_html
	args_copy = copy.deepcopy(args)
	args_copy.input_dir = args_copy.output_dir
	args_copy.command = args.cmd
	make_html.main(args_copy)
	
def chrombpnet_qc(args):

	if args.file_prefix:
		fpx = args.file_prefix+"_"
	else:
		fpx = ""
	
	def load_model_wrapper(model_hdf5):
		# read .h5 model
		from tensorflow.keras.utils import get_custom_objects
		from tensorflow.keras.models import load_model
		import tensorflow as tf
		import chrombpnet.training.utils.losses as losses
		custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
		get_custom_objects().update(custom_objects)    
		model=load_model(model_hdf5)
		model.summary()
		return model
    
	chrombpnet_md = load_model_wrapper(model_hdf5=args.chrombpnet_model)
	args.inputlen = int(chrombpnet_md.input_shape[1])
	args.outputlen = int(chrombpnet_md.output_shape[0][1])
	
	# make predictions with trained chrombpnet model
	import chrombpnet.training.predict as predict
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args.output_dir,"evaluation/{}chrombpnet".format(fpx))
	args_copy.model_h5 = args.chrombpnet_model
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/{}filtered.peaks.bed".format(fpx))
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
	args_copy.model_h5 = args.chrombpnet_model_nb
	args_copy.regions = args.nonpeaks
	args_copy.output_prefix = os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias".format(fpx))
	args_copy.motifs_to_pwm = os.path.join(args_copy.output_dir,"auxiliary/motif_to_pwm.tsv")
	args_copy.ylim = None
	marginal_footprinting.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias_footprints.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/{}chrombpnet_nobias_footprints.h5".format(fpx)))

	# get contributions scores with model
	args_copy = copy.deepcopy(args)
	import chrombpnet.evaluation.interpret.interpret as interpret
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/{}filtered.peaks.bed".format(fpx))
	peaks = pd.read_csv(os.path.join(args_copy.peaks),sep="\t",header=None)
	if peaks.shape[0] > 30000:
		sub_peaks = peaks.sample(30000, random_state=1234)
	else:
		sub_peaks = peaks
	sub_peaks.to_csv(os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx)),sep="\t", header=False, index=False)
	os.makedirs(os.path.join(args.output_dir,"auxiliary/interpret_subsample/"), exist_ok=False)

	#args_copy.profile_or_counts = ["counts", "profile"]
	args_copy.profile_or_counts = ["profile"]
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx))	
	args_copy.model_h5 = args.chrombpnet_model_nb
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias".format(fpx))
	args_copy.debug_chr = None
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
	# modisco-lite pipeline
	
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias.profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_profile/"),meme_file)
	os.system(modisco_command)
	#modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias.counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)))
	#os.system(modisco_command)
	#modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_counts/"),meme_file)
	#os.system(modisco_command)

	import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
	#convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_counts/motifs.html"),os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias_counts.pdf".format(fpx)))
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_profile/motifs.html"),os.path.join(args.output_dir,"evaluation/{}chrombpnet_nobias_profile.pdf".format(fpx)))
	
	import chrombpnet.helpers.generate_reports.make_html as make_html
	args_copy = copy.deepcopy(args)
	args_copy.input_dir = args_copy.output_dir
	args_copy.command = args.cmd
	make_html.main(args_copy)
	
def train_bias_pipeline(args):

	if args.file_prefix:
		fpx = args.file_prefix+"_"
	else:
		fpx = ""
		
	# Shift bam and convert to bigwig
	import chrombpnet.helpers.preprocessing.reads_to_bigwig as reads_to_bigwig	
	args.output_prefix = os.path.join(args.output_dir,"auxiliary/{}data".format(fpx))
	# args.plus_shift = None
	# args.minus_shift = None
	# reads_to_bigwig.main(args)
	
	# QC bigwig
	import chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig as build_pwm_from_bigwig	
	args.bigwig = os.path.join(args.output_dir,"auxiliary/{}data_unstranded.bw".format(fpx))
	args.output_prefix = os.path.join(args.output_dir,"evaluation/{}bw_shift_qc".format(fpx))
	folds = json.load(open(args.chr_fold_path))
	assert(len(folds["valid"]) > 0) # validation list of chromosomes is empty
	args.chr = folds["valid"][0]
	args.pwm_width=24
	build_pwm_from_bigwig.main(args)
	

	# fetch hyperparameters for training
	import chrombpnet.helpers.hyperparameters.find_bias_hyperparams as find_bias_hyperparams
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/{}".format(fpx))
	find_bias_hyperparams.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"auxiliary/{}bias_model_params.tsv".format(fpx)),os.path.join(args.output_dir,"logs/{}bias_model_params.tsv".format(fpx)))
	os.rename(os.path.join(args.output_dir,"auxiliary/{}bias_data_params.tsv".format(fpx)),os.path.join(args.output_dir,"logs/{}bias_data_params.tsv".format(fpx)))
		
	# get model architecture path and train chromBPNet model
	import chrombpnet.training.models.bpnet_model as bpnet_model
	import chrombpnet.training.train as train
	args_copy = copy.deepcopy(args)
	if args_copy.architecture_from_file is None:
		args_copy.architecture_from_file = 	bpnet_model.__file__
	args_copy.peaks = "None"
	args_copy.nonpeaks = os.path.join(args_copy.output_dir,"auxiliary/{}filtered.bias_nonpeaks.bed".format(fpx))
	args_copy.output_prefix = os.path.join(args_copy.output_dir,"models/{}bias".format(fpx))
	args_copy.params = os.path.join(args_copy.output_dir,"logs/{}bias_model_params.tsv".format(fpx))
	train.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"models/{}bias.args.json".format(fpx)),os.path.join(args.output_dir,"logs/{}bias.args.json".format(fpx)))
	os.rename(os.path.join(args.output_dir,"models/{}bias.log".format(fpx)),os.path.join(args.output_dir,"logs/{}bias.log".format(fpx)))
	os.rename(os.path.join(args.output_dir,"models/{}bias.log.batch".format(fpx)),os.path.join(args.output_dir,"logs/{}bias.log.batch".format(fpx)))
#	os.rename(os.path.join(args.output_dir,"models/{}bias.#".format(fpx)),os.path.join(args.output_dir,"logs/{}bias.params.json".format(fpx)))

	if args.cmd_bias == "train":
		import chrombpnet.helpers.generate_reports.make_html_bias as make_html_bias
		args_copy = copy.deepcopy(args)
		args_copy.input_dir = args_copy.output_dir
		args_copy.command = args_copy.cmd_bias
		make_html_bias.main(args_copy) 
		print("Finished training! Exiting!")
		return
		
	# make predictions with trained bias model 
	import chrombpnet.training.predict as predict
	args_copy = copy.deepcopy(args)
	args_copy.nonpeaks = os.path.join(args_copy.output_dir,"auxiliary/{}filtered.bias_nonpeaks.bed".format(fpx))
	args_copy.peaks = os.path.join(args_copy.output_dir,"auxiliary/{}filtered.bias_peaks.bed".format(fpx))
	args_copy.output_prefix = os.path.join(args_copy.output_dir,"evaluation/{}bias".format(fpx))
	args_copy.model_h5 = os.path.join(args.output_dir,"models/{}bias.h5".format(fpx))
	predict.main(args_copy)

	# get contributions scores with model
	import chrombpnet.evaluation.interpret.interpret as interpret
	args_copy.peaks = os.path.join(args_copy.output_dir,"auxiliary/{}filtered.bias_peaks.bed".format(fpx))
	peaks = pd.read_csv(os.path.join(args_copy.peaks),sep="\t",header=None)
	if peaks.shape[0] > 30000:
		sub_peaks = peaks.sample(30000, random_state=1234)
	else:
		sub_peaks = peaks
	sub_peaks.to_csv(os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx)),sep="\t", header=False, index=False)
	os.makedirs(os.path.join(args.output_dir,"auxiliary/interpret_subsample/"), exist_ok=False)

	args_copy = copy.deepcopy(args)
	args_copy.profile_or_counts = ["counts", "profile"]
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx))	
	args_copy.model_h5 = os.path.join(args.output_dir,"models/{}bias.h5".format(fpx))
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias".format(fpx))
	args_copy.debug_chr = None
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	# modisco-lite pipeline
	
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias.profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_profile/"),meme_file)
	os.system(modisco_command)
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias.counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_counts/"),meme_file)
	os.system(modisco_command)
	
	import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_counts/motifs.html"),os.path.join(args.output_dir,"evaluation/{}bias_counts.pdf".format(fpx)))
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_profile/motifs.html"),os.path.join(args.output_dir,"evaluation/{}bias_profile.pdf".format(fpx)))

	import chrombpnet.helpers.generate_reports.make_html_bias as make_html_bias
	args_copy = copy.deepcopy(args)
	args_copy.input_dir = args_copy.output_dir
	args_copy.command = args_copy.cmd_bias
	make_html_bias.main(args_copy)

def bias_model_qc(args):

	if args.file_prefix:
		fpx = args.file_prefix+"_"
	else:
		fpx = ""
	
	def load_model_wrapper(model_hdf5):
		# read .h5 model
		from tensorflow.keras.utils import get_custom_objects
		from tensorflow.keras.models import load_model
		import tensorflow as tf
		import chrombpnet.training.utils.losses as losses
		custom_objects={"multinomial_nll":losses.multinomial_nll, "tf": tf}    
		get_custom_objects().update(custom_objects)    
		model=load_model(model_hdf5)
		model.summary()
		return model
    
	bias_md = load_model_wrapper(model_hdf5=args.bias_model)
	args.inputlen = int(bias_md.input_shape[1])
	args.outputlen = int(bias_md.output_shape[0][1])
	
	# make predictions with trained bias model 
	import chrombpnet.training.predict as predict
	args_copy = copy.deepcopy(args)
	args_copy.nonpeaks = os.path.join(args_copy.output_dir,"auxiliary/{}filtered.bias_nonpeaks.bed".format(fpx))
	args_copy.peaks = os.path.join(args_copy.output_dir,"auxiliary/{}filtered.bias_peaks.bed".format(fpx))
	args_copy.output_prefix = os.path.join(args_copy.output_dir,"evaluation/{}bias".format(fpx))
	args_copy.model_h5 = args.bias_model
	predict.main(args_copy)

	# get contributions scores with model
	import chrombpnet.evaluation.interpret.interpret as interpret
	args_copy.peaks = os.path.join(args_copy.output_dir,"auxiliary/{}filtered.bias_peaks.bed".format(fpx))
	peaks = pd.read_csv(os.path.join(args_copy.peaks),sep="\t",header=None)
	if peaks.shape[0] > 30000:
		sub_peaks = peaks.sample(30000, random_state=1234)
	else:
		sub_peaks = peaks
	sub_peaks.to_csv(os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx)),sep="\t", header=False, index=False)
	os.makedirs(os.path.join(args.output_dir,"auxiliary/interpret_subsample/"), exist_ok=False)

	args_copy = copy.deepcopy(args)
	args_copy.profile_or_counts = ["counts", "profile"]
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx))	
	args_copy.model_h5 = args.bias_model
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias".format(fpx))
	args_copy.debug_chr = None
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
	# modisco-lite pipeline
	
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias.profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_profile/"),meme_file)
	os.system(modisco_command)
	modisco_command = "modisco motifs -i {} -n 50000 -o {} -w 500".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias.counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)))
	os.system(modisco_command)
	modisco_command = "modisco report -i {} -o {} -m {}".format(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_counts_scores.h5".format(fpx)),os.path.join(args.output_dir,"evaluation/modisco_counts/"),meme_file)
	os.system(modisco_command)
	
	import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_counts/motifs.html"),os.path.join(args.output_dir,"evaluation/{}bias_counts.pdf".format(fpx)))
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_profile/motifs.html"),os.path.join(args.output_dir,"evaluation/{}bias_profile.pdf".format(fpx)))

	import chrombpnet.helpers.generate_reports.make_html_bias as make_html_bias
	args_copy = copy.deepcopy(args)
	args_copy.input_dir = args_copy.output_dir
	args_copy.command = args_copy.cmd_bias
	make_html_bias.main(args_copy)

# --- internal posttrain stage helpers -------------------------------------
# Each `_pt_*_stage` runs one stage of the posttrain framework against an
# already-prepared output dir. They are called both by the end-to-end
# orchestrator (`posttrain_bias_correction_pipeline`) and by the standalone
# subcommand pipelines, so the same wiring isn't duplicated five times.


def _pt_setup_bigwig(args, fpx):
	"""
	Common bigwig prep + bigwig QC for any posttrain stage that needs the
	dataset's bigwig. Mutates `args`: sets `args.bigwig` to the prepared path
	and writes intermediate values to `args.output_prefix` along the way.

	If `args.input_bigwig` is set (via `--input-bigwig`), skip running
	`reads_to_bigwig` and use the supplied bigwig directly. The user is
	responsible for ensuring it's properly +4/-4 shifted.
	"""
	import chrombpnet.helpers.preprocessing.analysis.build_pwm_from_bigwig as build_pwm_from_bigwig
	if getattr(args, "input_bigwig", None):
		args.bigwig = args.input_bigwig
	else:
		import chrombpnet.helpers.preprocessing.reads_to_bigwig as reads_to_bigwig
		args.output_prefix = os.path.join(args.output_dir,"auxiliary/{}data".format(fpx))
		reads_to_bigwig.main(args)
		args.bigwig = os.path.join(args.output_dir,"auxiliary/{}data_unstranded.bw".format(fpx))

	args.output_prefix = os.path.join(args.output_dir,"evaluation/{}bw_shift_qc".format(fpx))
	folds = json.load(open(args.chr_fold_path))
	assert(len(folds["valid"]) > 0)
	args.chr = folds["valid"][0]
	args.pwm_width = 24
	build_pwm_from_bigwig.main(args)


def _pt_finetune_bias_stage(args, fpx):
	"""
	Filter background regions for the bias fine-tune step, then train the
	combined fine-tune model defined by `bpnet_finetune_combiner_model.py`.

	Pre: `args.bigwig` set; `args.pretrained_bias_model_paths` and
	`args.combiner_profile_kernel_size` populated by the CLI.
	Post: `models/<fpx>finetuned_bias.h5` produced.
	"""
	# `find_bias_hyperparams` reads filters/convgroups/n_dilation_layers/max_jitter
	# and writes them into bias_model_params.tsv. The combiner architecture
	# itself ignores those keys, but train.get_model_param_dict asserts they
	# exist, so we set them to standard bias defaults and force jitter=0.
	import chrombpnet.helpers.hyperparameters.find_bias_hyperparams as find_bias_hyperparams
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/{}".format(fpx))
	args_copy.filters = 128
	args_copy.convgroups = 1
	args_copy.n_dilation_layers = 4
	args_copy.max_jitter = 0
	find_bias_hyperparams.main(args_copy)

	os.rename(
		os.path.join(args.output_dir,"auxiliary/{}bias_model_params.tsv".format(fpx)),
		os.path.join(args.output_dir,"logs/{}bias_finetune_model_params.tsv".format(fpx)),
	)
	os.rename(
		os.path.join(args.output_dir,"auxiliary/{}bias_data_params.tsv".format(fpx)),
		os.path.join(args.output_dir,"logs/{}bias_finetune_data_params.tsv".format(fpx)),
	)
	finetune_params_path = os.path.join(args.output_dir,"logs/{}bias_finetune_model_params.tsv".format(fpx))
	with open(finetune_params_path, "a") as f:
		f.write("\t".join(["pretrained_bias_model_paths", ",".join(args.pretrained_bias_model_paths)]))
		f.write("\n")
		f.write("\t".join(["combiner_profile_kernel_size", str(args.combiner_profile_kernel_size)]))
		f.write("\n")

	import chrombpnet.training.models.bpnet_finetune_combiner_model as bpnet_finetune_combiner_model
	import chrombpnet.training.train as train
	args_copy = copy.deepcopy(args)
	if args_copy.architecture_from_file is None:
		args_copy.architecture_from_file = bpnet_finetune_combiner_model.__file__
	args_copy.peaks = "None"
	args_copy.nonpeaks = os.path.join(args.output_dir,"auxiliary/{}filtered.bias_nonpeaks.bed".format(fpx))
	args_copy.output_prefix = os.path.join(args.output_dir,"models/{}finetuned_bias".format(fpx))
	args_copy.params = finetune_params_path
	args_copy.max_jitter = 0
	train.main(args_copy)

	os.rename(
		os.path.join(args.output_dir,"models/{}finetuned_bias.args.json".format(fpx)),
		os.path.join(args.output_dir,"logs/{}finetuned_bias.args.json".format(fpx)),
	)
	os.rename(
		os.path.join(args.output_dir,"models/{}finetuned_bias.log".format(fpx)),
		os.path.join(args.output_dir,"logs/{}finetuned_bias.log".format(fpx)),
	)
	os.rename(
		os.path.join(args.output_dir,"models/{}finetuned_bias.log.batch".format(fpx)),
		os.path.join(args.output_dir,"logs/{}finetuned_bias.log.batch".format(fpx)),
	)


def _pt_filter_peaks_and_optionally_scale_bias(args, fpx, bias_model_path,
		bpnet_filters, bpnet_convgroups, bpnet_n_dilation_layers):
	"""
	Filter peaks/nonpeaks and write the BPNet-stage params TSV; if
	`bias_model_path != "None"`, also δ-scale the supplied bias model and
	move the scaled model to `models/`.

	Pre: `args.bigwig` set.
	Post:
	  - auxiliary/<fpx>filtered.peaks.bed, filtered.nonpeaks.bed
	  - logs/<fpx>bpnet_model_params.tsv, bpnet_data_params.tsv
	  - if bias_model_path != "None":
	      models/<fpx>bias_model_scaled.h5
	      bpnet_model_params.tsv's bias_model_path row points to that file
	Returns the path of the BPNet params TSV (for the BPNet train stage).
	"""
	import chrombpnet.helpers.hyperparameters.find_chrombpnet_hyperparams as find_chrombpnet_hyperparams
	args_copy = copy.deepcopy(args)
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/{}".format(fpx))
	args_copy.bias_model_path = bias_model_path
	args_copy.filters = bpnet_filters
	args_copy.convgroups = bpnet_convgroups
	args_copy.n_dilation_layers = bpnet_n_dilation_layers
	find_chrombpnet_hyperparams.main(args_copy)

	if bias_model_path != "None":
		os.rename(
			os.path.join(args.output_dir,"auxiliary/{}bias_model_scaled.h5".format(fpx)),
			os.path.join(args.output_dir,"models/{}bias_model_scaled.h5".format(fpx)),
		)
	os.rename(
		os.path.join(args.output_dir,"auxiliary/{}chrombpnet_model_params.tsv".format(fpx)),
		os.path.join(args.output_dir,"logs/{}bpnet_model_params.tsv".format(fpx)),
	)
	os.rename(
		os.path.join(args.output_dir,"auxiliary/{}chrombpnet_data_params.tsv".format(fpx)),
		os.path.join(args.output_dir,"logs/{}bpnet_data_params.tsv".format(fpx)),
	)

	bpnet_params_path = os.path.join(args.output_dir,"logs/{}bpnet_model_params.tsv".format(fpx))
	if bias_model_path != "None":
		params_text = open(bpnet_params_path).read().replace(
			os.path.join(args.output_dir,"auxiliary/{}bias_model_scaled.h5".format(fpx)),
			os.path.join(args.output_dir,"models/{}bias_model_scaled.h5".format(fpx)),
		)
		with open(bpnet_params_path, "w") as f:
			f.write(params_text)
	return bpnet_params_path


def _pt_train_bpnet_stage(args, fpx, params_path):
	"""
	Train the vanilla BPNet on the filtered peaks/nonpeaks using `params_path`.
	bpnet_model.py reads filters/convgroups/n_dil_layers/counts_loss_weight
	from the TSV; an `bias_model_path` row is harmlessly ignored.

	Pre: filtered.peaks.bed, filtered.nonpeaks.bed and params_path exist.
	Post: models/<fpx>bpnet.h5 produced.
	"""
	import chrombpnet.training.models.bpnet_model as bpnet_model
	import chrombpnet.training.train as train
	args_copy = copy.deepcopy(args)
	if args_copy.architecture_from_file is None:
		args_copy.architecture_from_file = bpnet_model.__file__
	args_copy.peaks = os.path.join(args.output_dir,"auxiliary/{}filtered.peaks.bed".format(fpx))
	args_copy.nonpeaks = os.path.join(args.output_dir,"auxiliary/{}filtered.nonpeaks.bed".format(fpx))
	args_copy.output_prefix = os.path.join(args.output_dir,"models/{}bpnet".format(fpx))
	args_copy.params = params_path
	train.main(args_copy)

	os.rename(
		os.path.join(args.output_dir,"models/{}bpnet.args.json".format(fpx)),
		os.path.join(args.output_dir,"logs/{}bpnet.args.json".format(fpx)),
	)
	os.rename(
		os.path.join(args.output_dir,"models/{}bpnet.log".format(fpx)),
		os.path.join(args.output_dir,"logs/{}bpnet.log".format(fpx)),
	)
	os.rename(
		os.path.join(args.output_dir,"models/{}bpnet.log.batch".format(fpx)),
		os.path.join(args.output_dir,"logs/{}bpnet.log.batch".format(fpx)),
	)


def _pt_subtract_stage(bpnet_h5, scaled_bias_h5, output_h5):
	import chrombpnet.training.models.chrombpnet_subtract_bias_model as chrombpnet_subtract_bias_model
	chrombpnet_subtract_bias_model.build_subtract_model(
		bpnet_h5=bpnet_h5, scaled_bias_h5=scaled_bias_h5, output_h5=output_h5,
	)


# --- public posttrain pipelines ------------------------------------------


def posttrain_bias_correction_pipeline(args):
	"""
	End-to-end post-training bias-correction.

	1. Fine-tune one or more pre-existing bias models on background regions
	   using the combiner architecture (`bpnet_finetune_combiner_model`).
	2. δ-scale the fine-tuned bias to the current dataset's read depth.
	3. Train a large vanilla BPNet on raw peaks (no bias attached).
	4. Build a Keras wrapper that subtracts the scaled bias from the BPNet
	   and save as `<prefix>chrombpnet_nobias.h5` so downstream tools
	   (`pred_bw`, `contribs_bw`, `footprints`, ...) work unchanged.
	"""
	fpx = args.file_prefix+"_" if args.file_prefix else ""
	print(args)

	_pt_setup_bigwig(args, fpx)
	_pt_finetune_bias_stage(args, fpx)

	finetuned_bias_h5 = os.path.join(args.output_dir,"models/{}finetuned_bias.h5".format(fpx))
	bpnet_params = _pt_filter_peaks_and_optionally_scale_bias(
		args, fpx,
		bias_model_path=finetuned_bias_h5,
		bpnet_filters=args.bpnet_filters,
		bpnet_convgroups=args.bpnet_convgroups,
		bpnet_n_dilation_layers=args.bpnet_n_dilation_layers,
	)
	_pt_train_bpnet_stage(args, fpx, params_path=bpnet_params)
	_pt_subtract_stage(
		bpnet_h5=os.path.join(args.output_dir,"models/{}bpnet.h5".format(fpx)),
		scaled_bias_h5=os.path.join(args.output_dir,"models/{}bias_model_scaled.h5".format(fpx)),
		output_h5=os.path.join(args.output_dir,"models/{}chrombpnet_nobias.h5".format(fpx)),
	)

	print("Post-training bias-correction pipeline complete.")
	print("Models produced under {}/models/:".format(args.output_dir))
	print("  {}finetuned_bias.h5".format(fpx))
	print("  {}bias_model_scaled.h5".format(fpx))
	print("  {}bpnet.h5".format(fpx))
	print("  {}chrombpnet_nobias.h5".format(fpx))


def posttrain_bpnet_pipeline(args):
	"""
	Standalone: train a vanilla BPNet on raw peaks (no bias attached). The
	model learns total signal (signal + bias). Pair its output with
	`posttrain bias-finetune` + `posttrain bias-scale` + `posttrain subtract`
	to get a bias-corrected `chrombpnet_nobias.h5` without retraining BPNet.
	"""
	fpx = args.file_prefix+"_" if args.file_prefix else ""
	print(args)

	_pt_setup_bigwig(args, fpx)
	# bias_model_path == "None" makes find_chrombpnet_hyperparams skip the
	# bias load+scale+save block and drop the bias_model_path TSV row.
	bpnet_params = _pt_filter_peaks_and_optionally_scale_bias(
		args, fpx,
		bias_model_path="None",
		bpnet_filters=args.bpnet_filters,
		bpnet_convgroups=args.bpnet_convgroups,
		bpnet_n_dilation_layers=args.bpnet_n_dilation_layers,
	)
	_pt_train_bpnet_stage(args, fpx, params_path=bpnet_params)
	print("BPNet trained on raw peaks. Output: {}/models/{}bpnet.h5".format(args.output_dir, fpx))


def posttrain_bias_finetune_pipeline(args):
	"""
	Standalone: fine-tune (and combine) one or more pre-existing bias models
	on background regions. Per pre-existing bias model, freeze every layer
	except the last dilated conv, the profile precrop conv, and the counts
	Dense; learn a Conv1D + Dense combiner on top of the heads.
	"""
	fpx = args.file_prefix+"_" if args.file_prefix else ""
	print(args)

	_pt_setup_bigwig(args, fpx)
	_pt_finetune_bias_stage(args, fpx)
	print("Combined bias model fine-tuned. Output: {}/models/{}finetuned_bias.h5".format(args.output_dir, fpx))


def posttrain_bias_scale_pipeline(args):
	"""
	Standalone: δ-shift the counts Dense bias of a bias model so its mean
	predicted log-counts on background regions matches mean(log1p(observed))
	on the current dataset. The shared helper also writes BPNet-stage params
	(harmless side outputs in `logs/`) and the filtered peaks/nonpeaks BEDs.
	"""
	fpx = args.file_prefix+"_" if args.file_prefix else ""
	print(args)

	_pt_setup_bigwig(args, fpx)
	# `args.bias_model_path` is the user-supplied bias to scale. The shared
	# helper's `bpnet_filters/convgroups/n_dilation_layers` end up only in the
	# side-output bpnet_model_params.tsv, which is unused for this command;
	# we pass the placeholder values that the bias-scale parser exposes.
	_pt_filter_peaks_and_optionally_scale_bias(
		args, fpx,
		bias_model_path=args.bias_model_path,
		bpnet_filters=args.filters,
		bpnet_convgroups=args.convgroups,
		bpnet_n_dilation_layers=args.n_dilation_layers,
	)
	print("Bias model scaled. Output: {}/models/{}bias_model_scaled.h5".format(args.output_dir, fpx))


def posttrain_subtract_pipeline(args):
	"""
	Standalone: build the Keras wrapper that produces bias-corrected
	predictions by subtracting a (scaled) bias model from a BPNet model.
	Saves a `*_nobias.h5` that downstream tools (`pred_bw`, `contribs_bw`,
	`footprints`, `modisco_motifs`) can consume directly.
	"""
	_pt_subtract_stage(args.bpnet_model_path, args.scaled_bias_model_path, args.output_path)
	print("Subtract wrapper saved: {}".format(args.output_path))
