import pandas as pd
import os
import json
import copy
from chrombpnet.data import DefaultDataFile, get_default_data_path
from chrombpnet.data import print_meme_motif_file
import numpy as np

def interpret_args(args_copy, args):
	# interpret reads args.seed/precision/batch_seqs; the pipeline's --seed and --precision are the training ones
	args_copy.seed = getattr(args, "shap_seed", 1234)
	args_copy.precision = getattr(args, "shap_precision", "auto")
	args_copy.batch_seqs = getattr(args, "shap_batch_seqs", None)
	return args_copy

def run_modisco(args, scores_h5, modisco_h5, report_dir, meme_file):
	# threads=None: the modisco subprocess gets the user's NUMBA_NUM_THREADS if set, else the CPUs allocated to the
	# job (chrombpnet.evaluation.modisco.run.default_threads)
	from chrombpnet.evaluation.modisco.run import modisco_motifs, modisco_report
	modisco_motifs(scores_h5, modisco_h5, max_seqlets=getattr(args, "modisco_max_seqlets", 50000),
		window=getattr(args, "modisco_window", 500), threads=None)
	modisco_report(modisco_h5, report_dir, str(meme_file), tomtom_lite=getattr(args, "tomtom_lite", False))

def run_modisco_heads(args, jobs, meme_file):
	# jobs: (scores_h5, modisco_h5, report_dir) per head, run one after the other like chrombpnet 1.x (each modisco
	# run uses all the CPUs of the job; running them at the same time would double the peak memory)
	for job in jobs:
		run_modisco(args, *job, meme_file)

def chrombpnet_train_pipeline(args):

	if args.file_prefix:
		fpx = args.file_prefix+"_"
	else:
		fpx = ""
		
	# Shift bam and convert to bigwig
	import chrombpnet.helpers.preprocessing.reads_to_bigwig as reads_to_bigwig	
	args.output_prefix = os.path.join(args.output_dir,"auxiliary/{}data".format(fpx))
	args.plus_shift = None
	args.minus_shift = None
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
	train.main(args_copy)
	
	# separating models from logs
	os.rename(os.path.join(args.output_dir,"models/{}chrombpnet.log".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet.log".format(fpx)))
	os.rename(os.path.join(args.output_dir,"models/{}chrombpnet.log.batch".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet.log.batch".format(fpx)))
	#os.rename(os.path.join(args.output_dir,"models/{}chrombpnet.params.json".format(fpx)),os.path.join(args.output_dir,"logs/{}chrombpnet.params.json".format(fpx)))
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
	interpret_subsample = getattr(args, "interpret_subsample", 30000)
	if peaks.shape[0] > interpret_subsample:
		sub_peaks = peaks.sample(interpret_subsample, random_state=1234)
	else:
		sub_peaks = peaks
	sub_peaks.to_csv(os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx)),sep="\t", header=False, index=False)
	os.makedirs(os.path.join(args.output_dir,"auxiliary/interpret_subsample/"), exist_ok=False)

	#args_copy.profile_or_counts = ["counts", "profile"]
	args_copy.profile_or_counts = ["profile"]	
	args_copy.regions = os.path.join(args_copy.output_dir,"auxiliary/{}30K_subsample_peaks.bed".format(fpx))	
	args_copy.model_h5 = os.path.join(args.output_dir,"models/{}chrombpnet_nobias.h5".format(fpx))
	args_copy.output_prefix = os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias".format(fpx))
	args_copy.debug_chr = None
	interpret_args(args_copy, args)
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
	# modisco-lite pipeline
	
	run_modisco(args, os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias.profile_scores.h5".format(fpx)),
		os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)),
		os.path.join(args.output_dir,"evaluation/modisco_profile/"), meme_file)
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
	
def chrombpnet_qc(args):

	if args.file_prefix:
		fpx = args.file_prefix+"_"
	else:
		fpx = ""
	
	from chrombpnet.training.utils.model_io import load_model_wrapper
	chrombpnet_md = load_model_wrapper(model_h5=args.chrombpnet_model)
	args.inputlen = int(chrombpnet_md.input_shape[1])
	args.outputlen = int(chrombpnet_md.output_shape[0][1])
	del chrombpnet_md
	
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
	interpret_subsample = getattr(args, "interpret_subsample", 30000)
	if peaks.shape[0] > interpret_subsample:
		sub_peaks = peaks.sample(interpret_subsample, random_state=1234)
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
	interpret_args(args_copy, args)
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
	# modisco-lite pipeline
	
	run_modisco(args, os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}chrombpnet_nobias.profile_scores.h5".format(fpx)),
		os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_profile_scores.h5".format(fpx)),
		os.path.join(args.output_dir,"evaluation/modisco_profile/"), meme_file)
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
	args.plus_shift = None
	args.minus_shift = None
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
	interpret_subsample = getattr(args, "interpret_subsample", 30000)
	if peaks.shape[0] > interpret_subsample:
		sub_peaks = peaks.sample(interpret_subsample, random_state=1234)
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
	interpret_args(args_copy, args)
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	# modisco-lite pipeline
	
	modisco_jobs = [(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias.{}_scores.h5".format(fpx, head)),
			os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_{}_scores.h5".format(fpx, head)),
			os.path.join(args.output_dir,"evaluation/modisco_{}/".format(head))) for head in ["profile", "counts"]]
	run_modisco_heads(args, modisco_jobs, meme_file)
	
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
	
	from chrombpnet.training.utils.model_io import load_model_wrapper
	bias_md = load_model_wrapper(model_h5=args.bias_model)
	args.inputlen = int(bias_md.input_shape[1])
	args.outputlen = int(bias_md.output_shape[0][1])
	del bias_md
	
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
	interpret_subsample = getattr(args, "interpret_subsample", 30000)
	if peaks.shape[0] > interpret_subsample:
		sub_peaks = peaks.sample(interpret_subsample, random_state=1234)
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
	interpret_args(args_copy, args)
	interpret.main(args_copy)
	
	import chrombpnet
	chrombpnet_src_dir = os.path.dirname(chrombpnet.__file__)
	meme_file=get_default_data_path(DefaultDataFile.motifs_meme)
	
	# modisco-lite pipeline
	
	modisco_jobs = [(os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}bias.{}_scores.h5".format(fpx, head)),
			os.path.join(args.output_dir,"auxiliary/interpret_subsample/{}modisco_results_{}_scores.h5".format(fpx, head)),
			os.path.join(args.output_dir,"evaluation/modisco_{}/".format(head))) for head in ["profile", "counts"]]
	run_modisco_heads(args, modisco_jobs, meme_file)
	
	import chrombpnet.evaluation.modisco.convert_html_to_pdf as convert_html_to_pdf
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_counts/motifs.html"),os.path.join(args.output_dir,"evaluation/{}bias_counts.pdf".format(fpx)))
	convert_html_to_pdf.main(os.path.join(args.output_dir,"evaluation/modisco_profile/motifs.html"),os.path.join(args.output_dir,"evaluation/{}bias_profile.pdf".format(fpx)))

	import chrombpnet.helpers.generate_reports.make_html_bias as make_html_bias
	args_copy = copy.deepcopy(args)
	args_copy.input_dir = args_copy.output_dir
	args_copy.command = args_copy.cmd_bias
	make_html_bias.main(args_copy)	
