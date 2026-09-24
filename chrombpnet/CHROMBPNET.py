import chrombpnet.parsers as parsers
import os
from chrombpnet.data import DefaultDataFile, get_default_data_path
from chrombpnet.data import print_meme_motif_file
import chrombpnet.pipelines as pipelines
import copy
import shutil
import subprocess
import pandas as pd
import logging
logging.getLogger('matplotlib.font_manager').disabled = True

MODEL_COMMANDS = ["pipeline", "train", "qc", "bias", "pred_bw", "contribs_bw", "footprints", "export"]

def runs_modisco(args):
	# the commands that end with modisco motifs + modisco report (train / bias train stop after training)
	if args.cmd == "bias":
		return args.cmd_bias in ("pipeline", "qc")
	return args.cmd in ("pipeline", "qc")

def check_tomtom(args):
	# modisco report needs MEME's tomtom unless --tomtom-lite; found the way evaluation/modisco/run.py finds it
	if not runs_modisco(args) or getattr(args, "tomtom_lite", False):
		return
	from chrombpnet.evaluation.modisco.run import modisco_env
	if shutil.which("tomtom", path=modisco_env()["PATH"]) is None:
		raise RuntimeError(
			"MEME's `tomtom` was not found on PATH; it is needed to match MoDISco motifs against {}. "
			"Install MEME (bioconda `meme`, included in chrombpnet's linux pixi environments) or use "
			"TOMTOM-lite instead (--tomtom-lite; reports p-values instead of q-values).".format(
				get_default_data_path(DefaultDataFile.motifs_meme)))

def check_model_runtime(args):
	# fail before any output directory is created
	check_tomtom(args)
	device = getattr(args, "device", None)
	if device == "cpu":
		# before jax is imported, so that JAX never creates a GPU client (a CUDA context on every visible GPU)
		os.environ["JAX_PLATFORMS"] = "cpu"
	if device in ("gpu", "cpu"):
		# 'auto' is left to the first model load (train.main logs the backend), so that the hours of CPU-only
		# preprocessing do not hold a context on every visible GPU
		from chrombpnet.training.runtime import assert_gpu_if_requested
		assert_gpu_if_requested(device)
	from chrombpnet.training.utils.model_io import require_jax_backend
	require_jax_backend()

def run_to_file(argv, output):
	# argv list, no shell: paths with spaces or shell metacharacters are passed to the command unchanged
	with open(output, "w") as out:
		subprocess.run(argv, stdout=out, check=True)

def bedtools_sort_merge(input_bed, output):
	# bedtools sort -i input_bed | bedtools merge -i stdin > output, failing if either command fails (like pipefail;
	# a failed sort would otherwise leave merge an empty stdin and an empty exclude.bed)
	with open(output, "w") as out, subprocess.Popen(["bedtools", "sort", "-i", input_bed], stdout=subprocess.PIPE) as sort:
		merge = subprocess.run(["bedtools", "merge", "-i", "stdin"], stdin=sort.stdout, stdout=out)
	# leaving the with block closed sort's stdout and waited for it
	for proc in (merge, sort):
		if proc.returncode != 0:
			raise subprocess.CalledProcessError(proc.returncode, proc.args)


# invoke pipeline modules based on command

def main():
	args = parsers.read_parser()

	if args.cmd in MODEL_COMMANDS:
		check_model_runtime(args)
	
	if args.cmd == "pipeline" or args.cmd == "train":
		os.makedirs(os.path.join(args.output_dir,"logs"), exist_ok=False)
		os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=False)
		os.makedirs(os.path.join(args.output_dir,"models"), exist_ok=False)
		os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=False)

		pipelines.chrombpnet_train_pipeline(args)
	
	elif args.cmd == "qc":
		os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=False)
		os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=False)
		
		pipelines.chrombpnet_qc(args)
		
	elif args.cmd == "bias":
		if args.cmd_bias == "pipeline" or args.cmd_bias == "train":
			os.makedirs(os.path.join(args.output_dir,"logs"), exist_ok=False)
			os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=False)
			os.makedirs(os.path.join(args.output_dir,"models"), exist_ok=False)
			os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=False)

			pipelines.train_bias_pipeline(args)
		
		elif args.cmd_bias == "qc":
			os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=False)
			os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=False)
			
			pipelines.bias_model_qc(args)
			
		else:
			print("Command not found")

			
	
	elif args.cmd == "pred_bw":
	
		assert (args.bias_model is None) + (args.chrombpnet_model is None) + (args.chrombpnet_model_nb is None) < 3, "No input model provided!"
		import chrombpnet.evaluation.make_bigwigs.predict_to_bigwig as predict_to_bigwig

		predict_to_bigwig.main(args)

	elif args.cmd == "contribs_bw":
	
		import chrombpnet.evaluation.interpret.interpret as interpret
		pipelines.interpret_args(args, args)
		interpret.main(args)
		import chrombpnet.evaluation.make_bigwigs.importance_hdf5_to_bigwig as importance_hdf5_to_bigwig
		if "counts" in  args.profile_or_counts:
			args_copy = copy.deepcopy(args)
			args_copy.hdf5 = args_copy.output_prefix + ".counts_scores.h5"
			args_copy.output_prefix = args.output_prefix + ".counts_scores"
			args_copy.regions =  args.output_prefix + ".interpreted_regions.bed"

			importance_hdf5_to_bigwig.main(args_copy)
		if "profile" in  args.profile_or_counts:
			args_copy = copy.deepcopy(args)
			args_copy.hdf5 = args_copy.output_prefix + ".profile_scores.h5"
			args_copy.output_prefix = args.output_prefix + ".profile_scores"
			args_copy.regions =  args.output_prefix + ".interpreted_regions.bed"
	
			importance_hdf5_to_bigwig.main(args_copy)
			
	elif args.cmd == "footprints":
	
		import chrombpnet.evaluation.marginal_footprints.marginal_footprinting as marginal_footprinting
		marginal_footprinting.main(args)

	elif args.cmd == "export":

		from chrombpnet.helpers.postprocessing.export_legacy_h5 import export_legacy_h5
		out = export_legacy_h5(args.model_h5, args.output, count_head=args.count_head)
		print("wrote {} ({}, count head: {})".format(out, args.format, args.count_head))

	elif args.cmd == "prep":
	
		if args.cmd_prep == "nonpeaks":

			assert(args.inputlen%2==0) # input length should be a multiple of 2
	
			os.makedirs(args.output_prefix+"_auxiliary/", exist_ok=False)
	
			from chrombpnet.helpers.make_gc_matched_negatives.get_genomewide_gc_buckets.get_genomewide_gc_bins import get_genomewide_gc
			get_genomewide_gc(args.genome,args.output_prefix+"_auxiliary/genomewide_gc.bed",args.inputlen, args.stride)
	
			# get gc content in peaks
			import chrombpnet.helpers.make_gc_matched_negatives.get_gc_content as get_gc_content
			args_copy = copy.deepcopy(args)
			args_copy.input_bed = args_copy.peaks
			args_copy.output_prefix = args.output_prefix+"_auxiliary/foreground.gc"
			get_gc_content.main(args_copy)
	
			# prepare candidate negatives
	
			exclude_bed = pd.read_csv(args.peaks, sep="\t", header=None)
			run_to_file(["bedtools", "slop", "-i", args.peaks, "-g", args.chrom_sizes, "-b", str(args.inputlen//2)],
					args.output_prefix+"_auxiliary/peaks_slop.bed")
			exclude_bed = pd.read_csv(args.output_prefix+"_auxiliary/peaks_slop.bed", sep="\t", header=None, usecols=[0,1,2])
	
			if args.blacklist_regions:
				run_to_file(["bedtools", "slop", "-i", args.blacklist_regions, "-g", args.chrom_sizes, "-b", str(args.inputlen//2)],
						args.output_prefix+"_auxiliary/blacklist_slop.bed")
										
				exclude_bed = pd.concat([exclude_bed,pd.read_csv(args.output_prefix+"_auxiliary/blacklist_slop.bed",sep="\t",header=None, usecols=[0,1,2])])

			exclude_bed.to_csv(args.output_prefix+"_auxiliary/exclude_unmerged.bed", sep="\t", header=False, index=False)
			bedtools_sort_merge(args.output_prefix+"_auxiliary/exclude_unmerged.bed", args.output_prefix+"_auxiliary/exclude.bed")

			run_to_file(["bedtools", "intersect", "-v", "-a", args.output_prefix+"_auxiliary/genomewide_gc.bed",
					"-b", args.output_prefix+"_auxiliary/exclude.bed"],
					args.output_prefix+"_auxiliary/candidates.bed")
													
			# get final negatives
			import chrombpnet.helpers.make_gc_matched_negatives.get_gc_matched_negatives as get_gc_matched_negatives
			args_copy = copy.deepcopy(args)
			args_copy.candidate_negatives = args.output_prefix+"_auxiliary/candidates.bed"
			args_copy.foreground_gc_bed = args.output_prefix+"_auxiliary/foreground.gc.bed"
			args_copy.output_prefix = 	args.output_prefix+"_auxiliary/negatives"
	
			get_gc_matched_negatives.main(args_copy)
	
			negatives = pd.read_csv(args.output_prefix+"_auxiliary/negatives.bed", sep="\t", header=None)
			negatives[3]="."
			negatives[4]="."
			negatives[5]="."
			negatives[6]="."
			negatives[7]="."
			negatives[8]="."
			negatives[9]=args.inputlen//2
			negatives.to_csv(args.output_prefix+"_negatives.bed", sep="\t", header=False, index=False)

		elif args.cmd_prep == "splits":
			import chrombpnet.helpers.make_chr_splits.splits as splits
			splits.main(args)
			
		else:
			print("Command not found")
		
	else:
		print("Command not found")


if __name__=="__main__":
	main()

    
        
