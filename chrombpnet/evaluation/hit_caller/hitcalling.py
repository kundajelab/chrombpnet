import os
import utils
import hitcaller
import hitcaller_background
import numpy as np
import pandas as pd
import subprocess
import argparse
import tqdm
import shutil
import qcutils
import matplotlib.pyplot as plt
import math

def fetch_arguments():
	parser=argparse.ArgumentParser(description="get hit calls with fimo or moods")
	parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
	parser.add_argument("-t", "--tomtom_anotation_path", type=str, required=False, default=None, help="If provided hits are renamed with Match_1 of tomtom annotations instead of modisco annotations.")
	parser.add_argument("-m", "--method", type=str, required=False, default="mean_norm", choices=['sum_norm', 'mean_norm', 'sum', 'mean'], help="Method to aggregate importance scores in hits ")
	parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Extracted peaks are centered at (2nd col) + summit (10th col).")
	parser.add_argument("-bgr", "--background_regions", type=str, required=True, help="10 column bed file of peaks. Extracted peaks are centered at (2nd col) + summit (10th col).")
	parser.add_argument("-w", "--width", type=int, required=False, default=1000, help="width of a peak around the summit. Hits are only called within these widths.")
	parser.add_argument("-cbw", "--contribs-bw", type=str, required=True, help="Bigwig file with contribution scores evaluated at the background regions bed file provided")
	parser.add_argument("-bgcbw", "--bg-contribs-bw", type=str, required=True, help="Bigwig file with contribution scores evaluated at the regions bed file provided")
	parser.add_argument("-mo", "--modisco-obj", type=str, required=True, help="Modisco object to extract pfms")
	parser.add_argument("--modisco-motifs-exclude", type=str, required=False, default=None, help="File path with modisco motif names (in format metacluster_pattern eg (0_0) per line.  If provied the motifs in file are excluded. If you have noisy motifs  (eg motif sequences of AAA or GGG) they degrade the quality of hitcaller so exclude them. You can also exlude hererodimers.")	
	parser.add_argument("-o", "--output-dir", type=str, required=True, help="Output directory to store results")
	parser.add_argument("-fdr","--contribs-fdr", type=float, required=False, default=0.05, help="FDR cut off to filter hits based on contribution scores")
	parser.add_argument("-pval","--init-scan-pval", type=str, required=False, default=0.01, help="pvalue cut off to use for mood or fimo for initial hitcalling")
	parser.add_argument("-hc","--hit-caller", type=str, required=False, default="fimo", choices=["moods", "fimo"], help="Use one of the hit callers")
	parser.add_argument("-th","--trim-threshold", type=float, required=False, default=0.3, help="trim threshold for motifs (Cut off anything less than x frac of max score)")
	parser.add_argument("--debug", action='store_true', default=False, help="If provided, all the intermediate files are not deleted after the run")
	args = parser.parse_args()
	return args

def get_thresholds_from_dinucs(args, pfms, cfms, suffix_in=""):
	
	#print(background_table.head())   (/0_0)
	suffix = "_bg"+suffix_in
	os.makedirs(os.path.join(args.output_dir,"auxiliary"+suffix), exist_ok=True) 

	# get hit calls in backgrounds with pfms
	suffix = "_bg"+suffix_in+"/pwms/"
	background_hit_table_pfms, pval_for_motif_pwm = hitcaller_background.get_hits(
		pfms, args.bg_contribs_bw,  hit_caller=args.hit_caller,
		temp_dir=os.path.join(args.output_dir,"auxiliary"+suffix), method=args.method, moods_pval_thresh=args.init_scan_pval
	)
	
	# get pval threshold
	pval_threshold_pwms = utils.compute_pval_threshold(background_hit_table_pfms, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir, suffix=suffix)

	# get hit calls in backgrounds with cwms
	suffix = "_bg"+suffix_in+"/cwms/"
	background_hit_table_cfms, pval_for_motif_cfms = hitcaller_background.get_hits(
		cfms, args.bg_contribs_bw, hit_caller=args.hit_caller,
		temp_dir=os.path.join(args.output_dir,"auxiliary"+suffix), method=args.method, moods_pval_thresh=args.init_scan_pval
	)
	
	# get pval threshold
	pval_threshold_cfms =  utils.compute_pval_threshold(background_hit_table_cfms, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir, suffix=suffix)
	
	return pval_threshold_pwms, pval_threshold_cfms, pval_for_motif_cfms, pval_for_motif_pwm
    

# def get_thresholds(args, pfms, cfms, suffix_in=""):
# 
# 	suffix = "_bg"
# 	background_bed_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix),"background_merged.bed")
# 
# 	if os.path.exists(background_bed_path):
# 		if os.stat(background_bed_path).st_size == 0:
# 			comm = ["bedtools"]
# 			comm += ["sort"]
# 			comm += ["-i"]
# 			comm += [args.background_regions]
# 			comm += ["|"]
# 			comm += ["bedtools"]
# 			comm += ["merge"]
# 			comm += ["-i"]
# 			comm += ["stdin"]
# 			f = open(background_bed_path, "w")
# 			proc = subprocess.Popen(" ".join(comm), shell=True, stdout=f)
# 			proc.wait()
# 	else:
# 		os.makedirs(os.path.join(args.output_dir,"auxiliary"+suffix), exist_ok=False) 
# 
# 		comm = ["bedtools"]
# 		comm += ["sort"]
# 		comm += ["-i"]
# 		comm += [args.background_regions]
# 		comm += ["|"]
# 		comm += ["bedtools"]
# 		comm += ["merge"]
# 		comm += ["-i"]
# 		comm += ["stdin"]
# 		f = open(background_bed_path, "w")
# 		proc = subprocess.Popen(" ".join(comm), shell=True, stdout=f)
# 		proc.wait()
# 			
# 	background_table = utils.import_peak_table_custom(background_bed_path)
# 		
# 	#print(background_table.head())   (/0_0)
# 	suffix = "_bg"+suffix_in
# 	os.makedirs(os.path.join(args.output_dir,"auxiliary"+suffix), exist_ok=True) 
# 
# 	# get hit calls in backgrounds with pfms
# 	suffix = "_bg"+suffix_in+"/pwms/"
# 	background_hit_table_pfms = hitcaller.get_hits(
# 		pfms, args.genome, background_bed_path, args.bg_contribs_bw, background_table, hit_caller=args.hit_caller,
# 		temp_dir=os.path.join(args.output_dir,"auxiliary"+suffix), method=args.method, moods_pval_thresh=args.init_scan_pval
# 	)
# 	
# 	# get pval threshold
# 	pval_threshold_pwms = utils.compute_pval_threshold(background_hit_table_pfms, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir, suffix=suffix)
# 
# 	# get hit calls in backgrounds with cwms
# 	suffix = "_bg"+suffix_in+"/cwms/"
# 	background_hit_table_cfms = hitcaller.get_hits(
# 		cfms, args.genome, background_bed_path, args.bg_contribs_bw, background_table, hit_caller=args.hit_caller,
# 		temp_dir=os.path.join(args.output_dir,"auxiliary"+suffix), method=args.method, moods_pval_thresh=args.init_scan_pval
# 	)
# 	
# 	# get pval threshold
# 	pval_threshold_cfms =  utils.compute_pval_threshold(background_hit_table_cfms, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir, suffix=suffix)
# 	
# 	return pval_threshold_pwms, pval_threshold_cfms
#     
def run_hit_caller(args, pfms, cfms, pval_threshold_pwms, pval_threshold_cfms, pval_for_motif_cfms, pval_for_motif_pwm, suffix_in=""):

	suffix = ""
	# center regions around the summit for peaks
	regions_bed_path=os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix),"summit_centered_regions.bed")
	peaks_bed_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix),"merged_peaks.bed")

	if os.path.exists(peaks_bed_path):
		if os.stat(peaks_bed_path).st_size == 0:
			regions = pd.read_csv(args.regions, sep="\t", header=None)
			regions[1] = regions[1] + regions[9] - args.width//2
			regions[2] = regions[1] +  args.width
			regions.to_csv(regions_bed_path, sep="\t", header=False, index=False)

			# merge peaks
			comm = ["bedtools"]
			comm += ["sort"]
			comm += ["-i"]
			comm += [os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix),"summit_centered_regions.bed")]
			comm += ["|"]
			comm += ["bedtools"]
			comm += ["merge"]
			comm += ["-i"]
			comm += ["stdin"]
			f = open(peaks_bed_path, "w")
			proc = subprocess.Popen(" ".join(comm), shell=True, stdout=f)
			proc.wait()
	else:
		os.makedirs(os.path.join(args.output_dir,"auxiliary"+suffix), exist_ok=False)

		regions = pd.read_csv(args.regions, sep="\t", header=None)
		regions[1] = regions[1] + regions[9] - args.width//2
		regions[2] = regions[1] +  args.width
		regions.to_csv(regions_bed_path, sep="\t", header=False, index=False)

		# merge peaks
		comm = ["bedtools"]
		comm += ["sort"]
		comm += ["-i"]
		comm += [os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix),"summit_centered_regions.bed")]
		comm += ["|"]
		comm += ["bedtools"]
		comm += ["merge"]
		comm += ["-i"]
		comm += ["stdin"]
		f = open(peaks_bed_path, "w")
		proc = subprocess.Popen(" ".join(comm), shell=True, stdout=f)
		proc.wait()
		
	peak_table = utils.import_peak_table_custom(peaks_bed_path)
	#print(peak_table.head())
	# run hit caller on the merged peaks

	suffix = suffix_in
	os.makedirs(os.path.join(args.output_dir,"auxiliary"+suffix), exist_ok=True) 

	suffix = suffix_in+"/pwms/"

	hit_table_pfms = hitcaller.get_hits(
			pfms, args.genome, peaks_bed_path, args.contribs_bw, peak_table, hit_caller=args.hit_caller,
			temp_dir=os.path.join(args.output_dir,"auxiliary"+suffix), method=args.method, moods_pval_thresh=pval_for_motif_pwm
		)

	hit_table_filtered_pfms = utils.use_threshold_for_filtering(hit_table_pfms, pval_threshold_pwms, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir, suffix=suffix)
	print("before filtering hits", hit_table_pfms.shape)
	print("after filtering hits", hit_table_filtered_pfms.shape)
	hit_table_filtered_pfms.to_csv(os.path.join(args.output_dir,"auxiliary"+suffix+"filtered.bed"), sep="\t", header=False, index=False)


	suffix = suffix_in+"/cwms/"
	hit_table_cfms = hitcaller.get_hits(
			cfms, args.genome, peaks_bed_path, args.contribs_bw, peak_table, hit_caller=args.hit_caller,
			temp_dir=os.path.join(args.output_dir,"auxiliary"+suffix), method=args.method, moods_pval_thresh=pval_for_motif_cfms
		)
	hit_table_cfms.to_csv(os.path.join(args.output_dir,"auxiliary"+suffix+"filtered.bed"), sep="\t", header=False, index=False)

	hit_table_filtered_cfms = utils.use_threshold_for_filtering(hit_table_cfms, pval_threshold_cfms, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir,  suffix=suffix)
	print("before filtering hits", hit_table_cfms.shape)
	print("after filtering hits", hit_table_filtered_cfms.shape)
	
	combine_filtered_hits = pd.concat([hit_table_filtered_pfms, hit_table_filtered_cfms], axis=0)

	#print(combine_filtered_hits.head())

	#hit_table_filtered = utils.filter_peak_hits_by_fdr(hit_table, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir, only_pos=only_pos,  suffix=suffix)

	suffix=suffix_in
	filtered_hits_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix), "moods_filtered.bed")
	combine_filtered_hits.to_csv(filtered_hits_path, sep="\t", header=False, index=False)
	

#	return new_hit_table

def resolve_overlaps(args, pfms, cfms, signs):

	suffix=""
	filtered_hits_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix), "moods_filtered.bed")

	# resolve overlaps
	hit_table = hitcaller.resolve_overlaps(os.path.join(args.output_dir,"auxiliary"+suffix),filtered_hits_path,args.contribs_bw,args.genome, args.tomtom_anotation_path, pfms, cfms, signs)
	hit_table.columns = ["chrom", "start", "end", "key", "strand", "match_score", "peak_index", "imp_score", "pvalue", "qvalue", "cluster", "cwm_activation_score"]

	filtered_hits_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix), "final.bed")
	hit_table.to_csv(filtered_hits_path, sep="\t", header=False, index=False)
	
	print("hits after resolving for overlaps", hit_table.shape)
	# annotate motifs with tomtom annotations
	if args.tomtom_anotation_path:
		tomtom = pd.read_csv(args.tomtom_anotation_path, sep="\t")
		label_dict = {}
		for index,row in tomtom.iterrows():
			keyd = str(row['Pattern']).replace("metacluster_","").replace("pattern_","").replace(".","_")
			label_dict[keyd] = keyd + "_" + str(row['Match_1'])
			hit_table['key'] = hit_table['key'].apply(lambda x: label_dict[x] if x in label_dict else x)


	new_hit_table = hit_table[["chrom", "start", "end", "key", "strand", "match_score", "imp_score", "cwm_activation_score", "pvalue", "qvalue"]]
	return new_hit_table


def qc_hit_caller(args):
	
	# intersect hit calls with peaks
	os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=True)
	os.makedirs(os.path.join(args.output_dir,"auxiliary_bg"), exist_ok=True)
	os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=True)

	eval_dir=os.path.join(args.output_dir,"evaluation")
	
	filtered_hits_path=os.path.join(os.path.join(args.output_dir,"auxiliary"),"hit_call_peak_indexed.bed")
	peak_bed_path=os.path.join(os.path.join(args.output_dir,"auxiliary"),"summit_centered_regions.bed")
	filter_hits_for_peaks = hitcaller.filter_hits_for_peaks(os.path.join(args.output_dir, "hit_calls.bed"), filtered_hits_path, peak_bed_path, 10)	
	
	hit_table_filtered = pd.read_csv(filtered_hits_path, sep="\t", header=None, names=["chrom", "start", "end", "key", "strand", "match_score", "imp_score", "cwm_score", "pvalue", "qvalue", "peak_index"])
	peak_table = utils.import_peak_table_custom(peak_bed_path)
	peak_hits = qcutils.get_peak_hits(peak_table, hit_table_filtered)
	motif_keys = list(set(hit_table_filtered["key"].values.tolist()))
	peak_hit_counts = qcutils.get_peak_motif_counts(peak_hits, motif_keys)
	
	# get proportion of hits with peaks
	motifs_per_peak = np.array([len(hits) for hits in peak_hits])
	print("Number of peaks with 0 motif hits: %d" % np.sum(motifs_per_peak == 0))
	
	
	quants = [0, 0.25, 0.50, 0.75, 0.99, 1]
	header = ["Quantile", "Number of hits/peak"]
	rows=[["%.1f%%" % (q * 100), "%d" % v]for q, v in zip(quants, np.quantile(motifs_per_peak, quants))]
	df = pd.DataFrame(rows,columns =header)
	df.to_csv(os.path.join(eval_dir,"hit_density_distribution.csv"), header=True, index=False)
	
	# hit density distribution
	fig, ax = plt.subplots(figsize=(10, 10))
	bins = np.concatenate([np.arange(np.max(motifs_per_peak) + 1), [np.inf]])
	ax.hist(motifs_per_peak, bins=bins, density=True, histtype="step", cumulative=True)
	ax.set_title("Cumulative distribution of number of motif hits per peak")
	ax.set_xlabel("Number of motifs k in peak")
	ax.set_ylabel("Proportion of peaks with at least k motifs")
	plt.savefig(os.path.join(eval_dir,"hit_density_distribution.png"))

	# hit versus motifs
	frac_peaks_with_motif = np.sum(peak_hit_counts > 0, axis=0) / len(peak_hit_counts)
	labels = np.array(motif_keys)
	sorted_inds = np.flip(np.argsort(frac_peaks_with_motif))
	frac_peaks_with_motif = frac_peaks_with_motif[sorted_inds]
	labels = labels[sorted_inds]
	
	fig, ax = plt.subplots(figsize=(15, 20))
	ax.bar(np.arange(len(labels)), frac_peaks_with_motif)
	ax.set_title("Proportion of peaks with each motif")
	ax.set_xticks(np.arange(len(labels)))
	ax.set_xticklabels(labels)
	plt.xticks(rotation=90)
	plt.savefig(os.path.join(eval_dir,"peaks_per_motif_hit.png"))

	# get cooccurence matrix 
	motif_cooccurrence_count_matrix = qcutils.get_motif_cooccurrence_count_matrix(peak_hit_counts)
	motif_cooccurrence_pval_matrix = qcutils.compute_cooccurrence_pvals(peak_hit_counts)
	
	per_motif_prob_cooccurence = motif_cooccurrence_count_matrix / np.sum(peak_hit_counts, axis=0)[:,None]
	
	qcutils.plot_peak_motif_indicator_heatmap(per_motif_prob_cooccurence, motif_keys, os.path.join(eval_dir,"per_motif_cooccurence_probability.png"))
	#store important objects for further analysis
	
	data = {"motifs_per_peak": motifs_per_peak,
	"peak_hit_counts": peak_hit_counts,
	"motif_keys": motif_keys,
	"motif_cooccurrence_count_matrix": motif_cooccurrence_count_matrix
	}
	
	np.save(os.path.join(eval_dir,'motif_stats.npy'),data)

	
if __name__=="__main__":

	args = fetch_arguments()
	exclude = open(args.modisco_motifs_exclude,"r").readlines()
	exclude = [x.strip() for x in exclude]
	
	if args.hit_caller == "moods" :
		out = shutil.which("moods-dna.py")
		assert out is not None, "Moods is not installed. Install moods and check  that command \"moods-dna.py -h\" returns a valid output without any error "
		
	if args.hit_caller == "fimo" :
		out = shutil.which("fimo")
		assert out is not None, "FIMO is not installed. Install FIMO and check  that command \"fimo -h\" returns a valid output without any error "

	# step 1 - prepare motif list for scanning
	init_pfms, init_cfms, sign_cfms, init_eq_cfms, init_eq_pfms = utils.import_tfmodisco_motifs(args.modisco_obj, exclude,  args.trim_threshold, only_pos=False)
	
	init_keys = list(set(list(init_pfms.keys()) + list(init_cfms.keys())))
	pfms = {}
	cfms = {}
	signs = {}
	eq_cfms = {}
	eq_pfms = {}
	
	motifs_keys = []
	
	print(args.modisco_motifs_exclude)
	print("inital motifs list used for hit calling: ", init_keys)
	if args.modisco_motifs_exclude is not None:
		exclude = open(args.modisco_motifs_exclude,"r").readlines()
		exclude = [x.strip() for x in exclude]
		print("provided list of motifs to exclude: ", exclude)
		motifs_keys = []
		for key in init_keys:
			if key in exclude:
				continue
			else:
				motifs_keys.append(key)
				pfms[key] = init_pfms[key]
				cfms[key] = init_cfms[key]
				signs[key] = sign_cfms[key]
				eq_cfms[key] = init_eq_cfms[key]
				eq_pfms[key] = init_eq_pfms[key]

	else:
		pfms = init_pfms
		cfms = init_cfms
		signs = sign_cfms
		motifs_keys = init_keys
		eq_cfms = init_eq_cfms
		eq_pfms = init_eq_pfms

		print("provided list of motifs to exclude: None")
	print("list of motifs to be used for hit calling: ", motifs_keys)
	
	# get thresholds for background
# 	
	# for each key get threshold and use it
# 	run_each_pwm = 0
# 	if run_each_pwm:
# 		for key in pfms:
# 			new_pfms = {}
# 			new_cfms = {}
# 			new_pfms[key] = pfms[key]
# 			new_cfms[key] = cfms[key]
# 	
# 			#pval_threshold_pwms,  pval_threshold_cfms = get_thresholds(args, new_pfms, new_cfms, "/"+key)
# 			pval_threshold_pwms,  pval_threshold_cfms = get_thresholds_from_dinucs(args, new_pfms, new_cfms, "/"+key)
# 
# 			print(pval_threshold_pwms,pval_threshold_cfms)
# 			# use computed thresholds on foreground filtering and resolve overlaps	
# 			run_hit_caller(args, new_pfms, new_cfms, pval_threshold_pwms, pval_threshold_cfms, "/"+key)
# 		
# 		frames = []
# 		for key in pfms:
# 			filtered_hits_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+"/"+key), "moods_filtered.bed")
# 			frames.append(pd.read_csv(filtered_hits_path, sep="\t", header=None))
# 			
# 		combine_filtered_hits = pd.concat(frames, axis=0)
# 		filtered_hits_full_path = os.path.join(os.path.join(args.output_dir,"auxiliary"), "moods_filtered.bed")
# 		combine_filtered_hits.to_csv(filtered_hits_full_path, sep="\t", header=False, index=False)
# 	else:
		#pval_threshold_pwms,  pval_threshold_cfms = get_thresholds(args, pfms, cfms, "")

	pval_threshold_pwms,  pval_threshold_cfms, pval_for_motif_cfms, pval_for_motif_pwm = get_thresholds_from_dinucs(args, pfms, cfms, "")

	print(pval_threshold_pwms,pval_threshold_cfms, pval_for_motif_cfms, pval_for_motif_pwm)
	# use computed thresholds on foreground filtering and resolve overlaps	
	run_hit_caller(args, pfms, cfms, pval_threshold_pwms, pval_threshold_cfms, pval_for_motif_cfms, pval_for_motif_pwm, "")


	# qc hit calls

	print("Resolve overlaps...")
	
	new_hit_table = resolve_overlaps(args, eq_pfms, eq_cfms, signs)
	new_hit_table.to_csv(os.path.join(args.output_dir, "hit_calls.bed"), sep="\t", header=False, index=False)

	qc_hit_caller(args)

	if not args.debug:
		shutil.rmtree(os.path.join(args.output_dir,"auxiliary"), ignore_errors=True)




