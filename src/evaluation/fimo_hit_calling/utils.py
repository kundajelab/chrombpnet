import h5py
import numpy as np
import pandas as pd
import tqdm
import pyBigWig
import os
import pomegranate
import sklearn.cluster
import scipy.cluster.hierarchy
import scipy.stats
import matplotlib.pyplot as plt
import subprocess


def import_peak_table_custom(peak_bed_paths):
	tables = []
	table = pd.read_csv(
		peak_bed_paths, sep="\t", header=None,  # Infer compression
		names=[
			"chrom", "peak_start", "peak_end"], usecols=[0,1,2]
	)
	return table

def import_shap_scores_from_bigwig( shap_scores_bw, peak_table):
    bw=pyBigWig.open(shap_scores_bw)
    scores=[]
    shape = []
    for i,r in peak_table.iterrows():
        val = np.nan_to_num(bw.values(r["chrom"],r["peak_start"],r["peak_end"]))
        scores.append(val)
        shape.append(len(val))

    return scores


def import_tfmodisco_motifs(tfm_results_path, exclude, trim_threshold, trim=True, only_pos=True):
	"""
	Imports the PFMs to into a dictionary, mapping `(x, y)` to the PFM,
	where `x` is the metacluster index and `y` is the pattern index.
	Arguments:
		`tfm_results_path`: path to HDF5 containing TF-MoDISco results
		`out_dir`: where to save motifs
		`trim`: if True, trim the motif flanks based on information content
		`only_pos`: if True, only return motifs with positive contributions
	Returns the dictionary of PFMs.
	""" 
    
	def trim_motif_new(cwm, motif, trim_threshold=0.3, max_motif_width=None):
		"""
		Given the PFM and motif (both L x 4 arrays) (the motif could be the
		PFM itself), trims `motif` by cutting off flanks of low information
		content in `pfm`. `min_ic` is the minimum required information
		content. If specified this trimmed motif will be extended on either
		side by `pad` bases.
		If no base passes the `min_ic` threshold, then no trimming is done.
		"""
		#print(cwm.shape)
		#print(motif.shape)


		score = np.sum(np.abs(cwm), axis=1)
		trim_thresh = np.max(score) * trim_threshold  # Cut off anything less than 30% of max score
		pass_inds = np.where(score >= trim_thresh)[0]
		trimmed = motif[np.min(pass_inds): np.max(pass_inds) + 1]

		#while trimmed.shape[0] <= 3:
		#	trim_threshold = trim_threshold - 0.1
		#	trim_thresh = np.max(score) * trim_threshold  # Cut off anything less than 30% of max score
		#	pass_inds = np.where(score >= trim_thresh)[0]
		#	trimmed = motif[np.min(pass_inds): np.max(pass_inds) + 1]
		
		#if not trimmed.size:
		#	return motif
		
		if max_motif_width is None:
			return trimmed
		else:
			assert(trimmed.shape[0] <= max_motif_width)
			
			right_extend_bp = 0
			left_extend_bp = 0
			#print(motif.shape)
			if trimmed.shape[0] == max_motif_width:
				return trimmed, left_extend_bp, right_extend_bp
				
			if trimmed.shape[0] < max_motif_width:
				add_bases = max_motif_width - trimmed.shape[0]
				#print(add_bases)
				right_idx = np.min(pass_inds)
				left_idx = np.max(pass_inds) + 1
				
				while add_bases > 0:
					if right_idx > 0:
						right_idx = right_idx - 1
						add_bases = add_bases - 1
						right_extend_bp += 1
						if add_bases == 0:
							break
					if left_idx < cwm.shape[0]:
						left_idx = left_idx + 1
						add_bases = add_bases - 1
						left_extend_bp += 1
					if right_idx == 0 and left_idx==cwm.shape[0]:
						break
				
				trimmed = motif[right_idx: left_idx]
				return trimmed, left_extend_bp, right_extend_bp
			
 
   
	def softmax(x, temp=100):
		norm_x = x - np.mean(x,axis=1, keepdims=True)
		return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

	pfms = {}
	cfms = {}
	sign_cfms = {}
	equal_width_cfms = {}
	equal_width_pfms = {}
	
	max_length = 0
	with h5py.File(tfm_results_path, "r") as f:
		metaclusters = f["metacluster_idx_to_submetacluster_results"]
		num_metaclusters = len(metaclusters.keys())
		for metacluster_i, metacluster_key in enumerate(metaclusters.keys()):
			metacluster = metaclusters[metacluster_key]
			if "patterns" not in metacluster["seqlets_to_patterns_result"]:
				continue
			patterns = metacluster["seqlets_to_patterns_result"]["patterns"]
			num_patterns = len(patterns["all_pattern_names"][:])
			for pattern_i, pattern_name in enumerate(patterns["all_pattern_names"][:]):
				#if pattern_i > 2:
				#	continue
				try:
					pattern_name = pattern_name.decode()
				except:
					pattern_name = pattern_name
				pattern = patterns[pattern_name]
				pfm = pattern["sequence"]["fwd"][:]
				cwm = pattern["task0_contrib_scores"]["fwd"][:]
				
				#if metacluster_i!=0 or pattern_i!=0:
				#	continue
					
				if "%d_%d" % (metacluster_i,pattern_i) in exclude:
					continue
					
				equal_width_pfms["%d_%d" % (metacluster_i,pattern_i)] = pfm

				pfm = trim_motif_new(cwm, pfm, trim_threshold)
				cfm = trim_motif_new(cwm, cwm, trim_threshold)
				pfms["%d_%d" % (metacluster_i,pattern_i)] = pfm
				equal_width_cfms["%d_%d" % (metacluster_i,pattern_i)] = cwm
				
				if np.sum(cwm)<0:
					equal_width_cfms["%d_%d" % (metacluster_i,pattern_i)] = cwm
					cfms["%d_%d" % (metacluster_i,pattern_i)] = softmax(cfm*-1)
					sign_cfms["%d_%d" % (metacluster_i,pattern_i)] = -1
				else:
					cfms["%d_%d" % (metacluster_i,pattern_i)] = softmax(cfm)
					sign_cfms["%d_%d" % (metacluster_i,pattern_i)] = +1
				
				if cfms["%d_%d" % (metacluster_i,pattern_i)].shape[0] > max_length:
					max_length = cfms["%d_%d" % (metacluster_i,pattern_i)].shape[0]
			
		print("max_motif_length", max_length)	
		for key in equal_width_cfms:
			#print(key)
				
			cwm = equal_width_cfms[key]	
			pfm = equal_width_pfms[key]	
			
			
			cfm, left_bp, right_bp = trim_motif_new(cwm, cwm, trim_threshold, max_length)
			pfm, left_bp_1, right_bp_1 = trim_motif_new(cwm, pfm, trim_threshold, max_length)
			
			assert(left_bp==left_bp_1)
			assert(right_bp==right_bp_1)
		
			equal_width_pfms[key] = pfm

			if np.sum(cwm)<0:
				equal_width_cfms[key] = [softmax(cfm*-1),left_bp,right_bp]
			else:
				equal_width_cfms[key] = [softmax(cfm),left_bp,right_bp]
				
	
			if np.sum(cwm)<0:
				equal_width_pfms[key] = [pfm,left_bp,right_bp]
			else:
				equal_width_pfms[key] = [pfm,left_bp,right_bp]
				
				
				#break					
				#print(cfms["%d_%d" % (metacluster_i,pattern_i)])
				#print(pfm)
						
# 				if only_pos and np.sum(cwm)<0:
# 					if trim:
# 						pfm = trim_motif_new(cwm, pfm, trim_threshold)
# 						cfm = trim_motif_new(cwm, cwm, trim_threshold)
# 						pfms["%d_%d" % (metacluster_i,pattern_i)] = pfm
# 						cfms["%d_%d" % (metacluster_i,pattern_i)] = cfm
# 				elif not only_pos  and np.sum(cwm)<0:
# 					if trim:
# 						pfm = trim_motif_new(cwm, pfm, trim_threshold)
# 						cfm = trim_motif_new(cwm, cwm, trim_threshold)
# 						pfms["%d_%d" % (metacluster_i,pattern_i)] = pfm
# 						cfms["%d_%d" % (metacluster_i,pattern_i)] = cfm
# 				else:
# 					continue
				# Check that the contribution scores are overall positive
				#if only_pos and np.sum(pfm<0) >0:
				#	continue	

	return pfms, cfms, sign_cfms, equal_width_cfms, equal_width_pfms
      
def compute_pval_threshold(hit_table, fdr_cutoff=0.05, output_dir=None, suffix="pos"):
	"""
	Filters the table of peak hits 
	by the importance score fraction by fitting a mixture model to the score
	distribution, taking the exponential component, and then fitting a
	percentile-tightened exponential distribution to this component.
	p-values are computed using this null, and then the FDR-cutoff is applied
	using Benjamini-Hochberg.
	Returns a reduced hit table of the same format. This will also generate
	plots for the score distribution and the FDR cutoffs.
	"""
	def exponential_pdf(x_values, lamb):
		return lamb * np.exp(-lamb * x_values)
	def exponential_cdf(x_values, lamb):
		return 1 - np.exp(-lamb * x_values) 
	
	def estimate_mode(x_values, bins=200, levels=1):
		"""
		Estimates the mode of the distribution using `levels`
		iterations of histograms.
		"""
		hist, edges = np.histogram(x_values, bins=bins)
		bin_mode = np.argmax(hist)
		left_edge, right_edge = edges[bin_mode], edges[bin_mode + 1]
		if levels <= 1:
			return (left_edge + right_edge) / 2
		else:
			return estimate_mode(
				x_values[(x_values >= left_edge) & (x_values < right_edge)],
				bins=bins,
				levels=(levels - 1)
			)
	
	def fit_tight_exponential_dist(x_values, mode=0, percentiles=np.arange(0.05, 1, 0.05)):
		"""
		Given an array of x-values and a set of percentiles of the distribution,
		computes the set of lambda values for an exponential distribution if the
		distribution were fit to each percentile of the x-values. Returns an array
		of lambda values parallel to `percentiles`. The exponential distribution
		is assumed to have the given mean/mode, and all data less than this mode
		is tossed out when doing this computation.
		"""
		assert np.min(percentiles) >= 0 and np.max(percentiles) <= 1
		x_values = x_values[x_values >= mode]
		per_x_vals = np.percentile(x_values, percentiles * 100)
		return -np.log(1 - percentiles) / (per_x_vals - mode)

	list_of_sampled_groups = []
	
	for name, group in hit_table.groupby('key'):
		print(group.shape)
		if group.shape[0] > 500000:
			sampled_group = group.sample(500000)
		else:
			sampled_group = group
	
	list_of_sampled_groups.append(sampled_group)
	output = pd.concat(list_of_sampled_groups).reset_index(drop=True)

	scores = output["imp_frac_score"].values

	scores_finite = scores[np.isfinite(scores)]
	
	from statsmodels.distributions.empirical_distribution import ECDF
	ecdf = ECDF(scores_finite)
	
	return ecdf
	

# 	mode = estimate_mode(scores_finite)
# 
# 	# Fit mixture of models to scores (mode-shifted)
# 	over_mode_scores = scores_finite[scores_finite >= mode] - mode
# 	mixed_model = pomegranate.GeneralMixtureModel.from_samples(
# 		[
# 			pomegranate.ExponentialDistribution,
# 			pomegranate.NormalDistribution,
# 			pomegranate.NormalDistribution
# 		],
# 		3, over_mode_scores[:, None]
# 	)
# 	mixed_model = mixed_model.fit(over_mode_scores)
# 	mixed_model_exp_dist = mixed_model.distributions[0]
# 
# 	# Obtain a distribution of scores that belong to the exponential distribution
# 	exp_scores = over_mode_scores[mixed_model.predict(over_mode_scores[:, None]) == 0]
# 
# 	# Fit a tight exponential distribution based on percentiles
# 	lamb = np.max(fit_tight_exponential_dist(exp_scores))
# 
# 	# Plot score distribution and fit
# 	
# 	import matplotlib
# 	
# 	matplotlib.rcParams['pdf.fonttype'] = 42
# 	matplotlib.rcParams['ps.fonttype'] = 42
# 
# 	fig, ax = plt.subplots(nrows=3, figsize=(20, 20))
# 
# 	x = np.linspace(np.min(scores_finite), np.max(scores_finite), 200)[1:]  # Skip first bucket (it's usually very large
# 	mix_dist_pdf = mixed_model.probability(x)
# 	mixed_model_exp_dist_pdf = mixed_model_exp_dist.probability(x)
# 
# 	perc_dist_pdf = exponential_pdf(x, lamb)
# 	perc_dist_cdf = exponential_cdf(x, lamb)
# 
# 	# Plot mixed model
# 	ax[0].hist(over_mode_scores + mode, bins=500, density=True, alpha=0.3)
# 	ax[0].axvline(mode)
# 	ax[0].plot(x + mode, mix_dist_pdf, label="Mixed model")
# 	ax[0].plot(x + mode, mixed_model_exp_dist_pdf, label="Exponential component")
# 	ax[0].set_ylabel("Density")
# 	ax[0].legend()
# 
# 	# Plot fitted PDF
# 	ax[1].hist(exp_scores, bins=500, density=True, alpha=0.3)
# 	ax[1].plot(x + mode, perc_dist_pdf, label="Percentile-fitted")
# 	ax[1].set_ylabel("Density")
# 
# 	# Plot fitted CDF
# 	ax[2].hist(exp_scores, bins=500, density=True, alpha=1, cumulative=True, histtype="step")
# 	ax[2].plot(x + mode, perc_dist_cdf, label="Percentile-fitted")
# 	ax[2].set_xlabel("Contribution scores aggregated in hits")
# 	ax[2].set_ylabel("CDF")
# 	#plt.savefig(os.path.join(output_dir,"score_distribution_and_fit.pdf"), dpi=300, transparent=True)
# 	plt.savefig(os.path.join(os.path.join(output_dir, "auxiliary"+suffix),"score_distribution_and_fit.pdf"), dpi=300, transparent=True)
# 
# 	# Compute p-values
# 	score_range = np.linspace(np.min(scores_finite), np.max(scores_finite), 1000000)
# 	inverse_cdf = 1 - exponential_cdf(score_range, lamb)
# 	assignments = np.digitize(scores - mode, score_range, right=True)
# 	assignments[~np.isfinite(scores)] = 0  # If score was NaN, give it a p-value of ~1
# 	pvals = inverse_cdf[assignments]
# 	pvals_sorted = np.sort(pvals)
# 
# 	# Plot FDR cut-offs of various levels
# 	fdr_levels = [0.05, 0.1, 0.2, 0.3]
# 	pval_threshes = []
# 	fig, ax = plt.subplots(figsize=(20, 8))
# 	ranks = np.arange(1, len(pvals_sorted) + 1)
# 	ax.plot(ranks, pvals_sorted, color="black", label="p-values")
# 	for fdr in fdr_levels:
# 		bh_crit_vals = (ranks / len(ranks)) * fdr
# 		ax.plot(ranks, bh_crit_vals, label=("Crit values (FDR = %.2f)" % fdr))
# 		inds = np.where(pvals_sorted <= bh_crit_vals)[0]
# 		if not len(inds):
# 			pval_threshes.append(-1)
# 		else:
# 			pval_threshes.append(pvals_sorted[np.max(inds)])
# 	ax.set_title("Step-up p-values and FDR corrective critical values")
# 	plt.legend()
# 	#plt.show()
# 	plt.savefig(os.path.join(os.path.join(output_dir, "auxiliary"+suffix),"pvals_fdr_corrected.pdf"), dpi=300, transparent=True)
# 
# 	# Show table of number of hits at each FDR level
# 	header = ["FDR level","Number of hits kept", "% hits kept"]
# 	rows = []
# 	rows.append(["None", len(pvals) ,"100%"])
# 	for i, fdr in enumerate(fdr_levels):
# 		num_kept = np.sum(pvals <= pval_threshes[i])
# 		frac_kept = num_kept / len(pvals)
# 		rows.append(["{0:.2f}".format(fdr), "{0}".format(int(num_kept)), "{0:.2f}".format(frac_kept * 100)])
# 	df = pd.DataFrame(rows,columns =header)
# 	df.to_csv(os.path.join(os.path.join(output_dir, "auxiliary"+suffix),"fdr_versus_hits.csv"), header=True, index=False)
# 	
# 	# Perform filtering
# 	bh_crit_vals = fdr_cutoff * (ranks / len(ranks))
# 	inds = np.where(pvals_sorted <= bh_crit_vals)[0]
# 	if not len(inds):
# 		pval_thresh = -1
# 	else:
# 		pval_thresh = pvals_sorted[np.max(inds)]
# 	return (pval_thresh, lamb, mode)
#     
def use_threshold_for_filtering(hit_table, ecdf, fdr_cutoff, output_dir=None, suffix="pos"):
	"""
	Filters the table of peak hits 
	by the importance score fraction by fitting a mixture model to the score
	distribution, taking the exponential component, and then fitting a
	percentile-tightened exponential distribution to this component.
	p-values are computed using this null, and then the FDR-cutoff is applied
	using Benjamini-Hochberg.
	Returns a reduced hit table of the same format. This will also generate
	plots for the score distribution and the FDR cutoffs.
	"""
	
# 	thresh = pval_thresh[0]
# 	lamb = pval_thresh[1]
	#mode = pval_thresh[2]
# 	
# 	def exponential_pdf(x_values, lamb):
# 		return lamb * np.exp(-lamb * x_values)
# 	def exponential_cdf(x_values, lamb):
# 		return 1 - np.exp(-lamb * x_values) 
# 	
# 	def estimate_mode(x_values, bins=200, levels=1):
# 		"""
# 		Estimates the mode of the distribution using `levels`
# 		iterations of histograms.
# 		"""
# 		hist, edges = np.histogram(x_values, bins=bins)
# 		bin_mode = np.argmax(hist)
# 		left_edge, right_edge = edges[bin_mode], edges[bin_mode + 1]
# 		if levels <= 1:
# 			return (left_edge + right_edge) / 2
# 		else:
# 			return estimate_mode(
# 				x_values[(x_values >= left_edge) & (x_values < right_edge)],
# 				bins=bins,
# 				levels=(levels - 1)
# 			)
# 	
# 	def fit_tight_exponential_dist(x_values, mode=0, percentiles=np.arange(0.05, 1, 0.05)):
# 		"""
# 		Given an array of x-values and a set of percentiles of the distribution,
# 		computes the set of lambda values for an exponential distribution if the
# 		distribution were fit to each percentile of the x-values. Returns an array
# 		of lambda values parallel to `percentiles`. The exponential distribution
# 		is assumed to have the given mean/mode, and all data less than this mode
# 		is tossed out when doing this computation.
# 		"""
# 		assert np.min(percentiles) >= 0 and np.max(percentiles) <= 1
# 		x_values = x_values[x_values >= mode]
# 		per_x_vals = np.percentile(x_values, percentiles * 100)
# 		return -np.log(1 - percentiles) / (per_x_vals - mode)
# 
# 
	scores = hit_table["imp_frac_score"].values 
	scores_finite = scores[np.isfinite(scores)]
	assert(scores_finite.shape==scores.shape) # some of the scores are not finite
	pvals = 1-ecdf(scores_finite)
	
	pvals_sorted = np.sort(pvals)
	
	# Plot FDR cut-offs of various levels
	fdr_levels = [0.05, 0.1, 0.2, 0.3]
	pval_threshes = []
	fig, ax = plt.subplots(figsize=(20, 8))
	ranks = np.arange(1, len(pvals_sorted) + 1)
	ax.plot(ranks, pvals_sorted, color="black", label="p-values")
	for fdr in fdr_levels:
		bh_crit_vals = (ranks / len(ranks)) * fdr
		ax.plot(ranks, bh_crit_vals, label=("Crit values (FDR = %.2f)" % fdr))
		inds = np.where(pvals_sorted <= bh_crit_vals)[0]
		if not len(inds):
			pval_threshes.append(-1)
		else:
			pval_threshes.append(pvals_sorted[np.max(inds)])
	ax.set_title("Step-up p-values and FDR corrective critical values")
	plt.legend()
	#plt.show()
	plt.savefig(os.path.join(os.path.join(output_dir, "auxiliary"+suffix),"pvals_fdr_corrected.pdf"), dpi=300, transparent=True)

	# Show table of number of hits at each FDR level
	header = ["FDR level","Number of hits kept", "% hits kept"]
	rows = []
	rows.append(["None", len(pvals) ,"100%"])
	for i, fdr in enumerate(fdr_levels):
		num_kept = np.sum(pvals <= pval_threshes[i])
		frac_kept = num_kept / len(pvals)
		rows.append(["{0:.2f}".format(fdr), "{0}".format(int(num_kept)), "{0:.2f}".format(frac_kept * 100)])
	df = pd.DataFrame(rows,columns =header)
	print(df)
	df.to_csv(os.path.join(os.path.join(output_dir, "auxiliary"+suffix),"fdr_versus_hits.csv"), header=True, index=False)
	
	# Perform filtering
	bh_crit_vals = fdr_cutoff * (ranks / len(ranks))
	inds = np.where(pvals_sorted <= bh_crit_vals)[0]
	if not len(inds):
		pval_thresh = -1
	else:
		pval_thresh = pvals_sorted[np.max(inds)]
	
	
	from statsmodels.stats.multitest import fdrcorrection
	rejected,pvalue_corrected =  fdrcorrection(pvals,alpha=fdr_cutoff)
	#print(np.sum(rejected))
	#print(np.sum(pvals <= pval_thresh))
	
	hit_table["pvalue"] = pvals
	hit_table["qvalue"] = pvalue_corrected
	
	return hit_table.iloc[pvals <= pval_thresh]				
# 	mode = estimate_mode(scores_finite)
# 	
# 	# Compute p-values
# 	score_range = np.linspace(np.min(scores_finite), np.max(scores_finite), 1000000)
# 	
# 	inverse_cdf = 1 - exponential_cdf(score_range, lamb) # lambda from the null model
# 	
# 	assignments = np.digitize(scores - mode, score_range, right=True) # center around null models mode
# 	
# 	assignments[~np.isfinite(scores)] = 0  # If score was NaN, give it a p-value of ~1
# 	pvals = inverse_cdf[assignments]


	return hit_table.iloc[pvals <= thresh]
	
