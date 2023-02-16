import os
import utils
import hitcaller
import numpy as np
import pandas as pd
import subprocess
import argparse
import tqdm
import shutil
import qcutils
import matplotlib.pyplot as plt

def fetch_arguments():
	parser=argparse.ArgumentParser(description="get hit calls with fimo or moods")
	parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
	parser.add_argument("-t", "--tomtom_anotation_path", type=str, required=False, default=None, help="If provided hits are renamed with Match_1 of tomtom annotations instead of modisco annotations.")
	parser.add_argument("-m", "--method", type=str, required=False, default="mean_norm", choices=['sum_norm', 'mean_norm', 'sum', 'mean'], help="Method to aggregate importance scores in hits ")
	parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. Extracted peaks are centered at (2nd col) + summit (10th col).")
	parser.add_argument("-w", "--width", type=int, required=False, default=1000, help="width of a peak around the summit. Hits are only called within these widths.")
	parser.add_argument("-cbw", "--contribs-bw", type=str, required=True, help="Bigwig file with contribution scores evaluated at the regions bed file provided")
	parser.add_argument("-mo", "--modisco-obj", type=str, required=True, help="Modisco object to extract pfms")
	parser.add_argument("--modisco-motifs-exclude", type=str, required=False, default=None, help="File path with modisco motif names (in format metacluster_pattern eg (0_0) per line.  If provied the motifs in file are excluded. If you have noisy motifs  (eg motif sequences of AAA or GGG) they degrade the quality of hitcaller so exclude them. You can also exlude hererodimers.")	
	parser.add_argument("-o", "--output-dir", type=str, required=True, help="Output directory to store results")
	parser.add_argument("-fdr","--contribs-fdr", type=float, required=False, default=0.1, help="FDR cut off to filter hits based on contribution scores")
	parser.add_argument("-pval","--init-scan-pval", type=str, required=False, default=0.001, help="pvalue cut off to use for mood or fimo for initial hitcalling")
	parser.add_argument("-hc","--hit-caller", type=str, required=False, default="fimo", choices=["moods", "fimo"], help="Use one of the hit callers")
	parser.add_argument("-th","--trim-threshold", type=float, required=False, default=0.3, help="trim threshold for motifs (Cut off anything less than x frac of max score)")
	parser.add_argument("--debug", action='store_true', default=False, help="If provided, all the intermediate files are not deleted after the run")
	args = parser.parse_args()
	return args


def run_hit_caller(args, init_pfms, init_cfms, only_pos):

	suffix = ""
	
	init_keys = list(init_pfms.keys())
	pfms = {}
	cfms = {}
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

	else:
		pfms = init_pfms
		cfms = init_cfms
		motifs_keys = init_keys
		print("provided list of motifs to exclude: None")
	
	print("list of motifs to be used for hit calling: ", motifs_keys)
		
	
	os.makedirs(os.path.join(args.output_dir,"auxiliary"+suffix), exist_ok=False)
	regions = pd.read_csv(args.regions, sep="\t", header=None)
	regions[1] = regions[1] + regions[9] - args.width//2
	regions[2] = regions[1] +  args.width
	regions.to_csv(os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix),"summit_centered_regions.bed"), sep="\t", header=False, index=False)
	peaks_bed_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix),"merged_peaks.bed")

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
	print(peak_table.head())
	
	if args.hit_caller == "moods" :
		out = shutil.which("moods-dna.py")
		assert out is not None, "Moods is not installed. Install moods and check  that command \"moods-dna.py -h\" returns a valid output without any error "

	if args.hit_caller == "fimo" :
		out = shutil.which("fimo")
		assert out is not None, "FIMO is not installed. Install FIMO and check  that command \"fimo -h\" returns a valid output without any error "

	hit_table = hitcaller.get_hits(
			pfms, args.genome, peaks_bed_path, args.contribs_bw, peak_table, hit_caller=args.hit_caller,
			temp_dir=os.path.join(args.output_dir,"auxiliary"+suffix), method=args.method, moods_pval_thresh=args.init_scan_pval
		)

	print(hit_table.head())
	print("Filter hits based on FDR...")

	# Filter motif hit table by p-value using FDR estimation
	hit_table_filtered = utils.filter_peak_hits_by_fdr(hit_table, fdr_cutoff=args.contribs_fdr, output_dir=args.output_dir, only_pos=only_pos,  suffix=suffix)

	filtered_hits_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix), "moods_filtered.bed")
	hit_table_filtered.to_csv(filtered_hits_path, sep="\t", header=False, index=False)
	
	print("Resolve overlaps...")

	# resolve overlaps
	hit_table = hitcaller.resolve_overlaps(os.path.join(args.output_dir,"auxiliary"+suffix),filtered_hits_path,args.contribs_bw,args.genome, pfms)
	hit_table.columns = ["chrom", "start", "end", "key", "strand", "match_score", "peak_index", "imp_score", "cluster", "cwm_activation_score"]

	# annotate motifs with tomtom annotations
	if args.tomtom_anotation_path:
		tomtom = pd.read_csv(args.tomtom_anotation_path, sep="\t")
		label_dict = {}
		for index,row in tomtom.iterrows():
			keyd = str(row['Pattern']).replace("metacluster_","").replace("pattern_","").replace(".","_")
			label_dict[keyd] = keyd + "_" + str(row['Match_1'])
			hit_table['key'] = hit_table['key'].apply(lambda x: label_dict[x] if x in label_dict else x)

	filtered_hits_path = os.path.join(os.path.join(args.output_dir,"auxiliary"+suffix), "final.bed")
	hit_table.to_csv(filtered_hits_path, sep="\t", header=False, index=False)

	new_hit_table = hit_table[["chrom", "start", "end", "key", "strand", "imp_score", "cwm_activation_score"]]

	return new_hit_table


def qc_hit_caller(args):
	
	# intersect hit calls with peaks
	os.makedirs(os.path.join(args.output_dir,"evaluation"), exist_ok=True)
	os.makedirs(os.path.join(args.output_dir,"auxiliary"), exist_ok=True)

	eval_dir=os.path.join(args.output_dir,"evaluation")
	
	filtered_hits_path=os.path.join(os.path.join(args.output_dir,"auxiliary"),"hit_call_peak_indexed.bed")
	peak_bed_path=os.path.join(os.path.join(args.output_dir,"auxiliary"),"summit_centered_regions.bed")
	filter_hits_for_peaks = hitcaller.filter_hits_for_peaks(os.path.join(args.output_dir, "hit_calls.bed"), filtered_hits_path, peak_bed_path, 10)	
	
	hit_table_filtered = pd.read_csv(filtered_hits_path, sep="\t", header=None, names=["chrom", "start", "end", "key", "strand", "imp_score", "cwm_score", "peak_index"])
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
	
	
	init_pfms, init_cfms = utils.import_tfmodisco_motifs(args.modisco_obj, args.trim_threshold, only_pos=False)

	#run_hit_caller(args)
	# run hit caller for positive contribution motifs
	#new_hit_table_pos = run_hit_caller(args, init_pfms_pos, init_cfms_pos, only_pos=True)
	# run hit caller for negative contribution motifs
	new_hit_table = run_hit_caller(args, init_pfms, init_cfms, only_pos=False)
	
	#if new_hit_table_neg is not None:
	#	new_hit_table_neg.to_csv(os.path.join(args.output_dir, "hit_calls_neg.bed"), sep="\t", header=False, index=False)

	#if new_hit_table_pos is not None:
	#	new_hit_table_pos.to_csv(os.path.join(args.output_dir, "hit_calls_pos.bed"), sep="\t", header=False, index=False)

	#print(new_hit_table_pos.shape)
	
	#new_hit_table_neg = pd.read_csv(os.path.join(args.output_dir, "hit_calls_neg.bed"), sep="\t", header=None)
	print(new_hit_table.shape)

	#new_hit_table = pd.concat([new_hit_table_neg, new_hit_table_pos], axis=0)

	new_hit_table.to_csv(os.path.join(args.output_dir, "hit_calls.bed"), sep="\t", header=False, index=False)

	qc_hit_caller(args)

	if not args.debug:
		shutil.rmtree(os.path.join(args.output_dir,"auxiliary"), ignore_errors=True)




