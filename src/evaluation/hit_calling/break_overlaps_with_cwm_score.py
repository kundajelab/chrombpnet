
import pandas as pd
from tqdm import tqdm
import argparse
import subprocess
import os
import pyfaidx
import pyBigWig
import h5py
import numpy as np
import one_hot
import scipy.stats
from scipy import signal

parser = argparse.ArgumentParser(description="get sequence contribution scores for the model")
parser.add_argument("-o", "--output_prefix", type=str, required=True, help="Output prefix")
parser.add_argument("-tfm", "--tfm_path", type=str, required=True, help="TFM path")
parser.add_argument("-bw", "--bigwig", type=str, required=True, help="bigwig path")
parser.add_argument("-g", "--genome", type=str, required=True, help="genome path")
args = parser.parse_args()



#inpath = os.path.join(args.output_prefix,"moods_filtered_scored_thresholded_normed.bed")

#comm = ["bedtools", "sort", "-i"]
#comm += [inpath]
#comm += ["|"]
#comm += ["bedtools", "cluster", "-i", "stdin"]
#f = open(os.path.join(args.output_prefix,"moods_filtered_scored_thresholded_clustered_normed.bed"), "w")
#proc = subprocess.Popen(" ".join(comm), shell=True, stdout=f)
#proc.wait()
#f.close()

def trim_motif_new(cwm, motif, trim_threshold=0.3):
    """
    Given the PFM and motif (both L x 4 arrays) (the motif could be the
    PFM itself), trims `motif` by cutting off flanks of low information
    content in `pfm`. `min_ic` is the minimum required information
    content. If specified this trimmed motif will be extended on either
    side by `pad` bases.
    If no base passes the `min_ic` threshold, then no trimming is done.
    """
    
    score = np.sum(np.abs(cwm), axis=1)
    trim_thresh = np.max(score) * trim_threshold  # Cut off anything less than 30% of max score
    pass_inds = np.where(score >= trim_thresh)[0]
    trimmed = motif[np.min(pass_inds): np.max(pass_inds) + 1]
 
    if not trimmed.size:
        return motif
    
    return trimmed

def import_tfmodisco_motifs(tfm_results_path, trim=True, only_pos=True):
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
    pfms = {}
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
                pattern_name = pattern_name.decode()
                pattern = patterns[pattern_name]
                pfm = pattern["sequence"]["fwd"][:]
                cwm = pattern["task0_contrib_scores"]["fwd"][:]
                
                # Check that the contribution scores are overall positive
                if only_pos and np.sum(cwm) < 0:
                    continue
                    
                if trim:
                    pfm = trim_motif_new(cwm, pfm)
                    
                pfms["%d_%d" % (metacluster_i,pattern_i)] = pfm
    return pfms

pfms = import_tfmodisco_motifs(args.tfm_path)


#data = pd.read_csv(os.path.join(args.output_prefix,"moods_filtered_scored_thresholded_clustered_subset.bed"), sep="\t", header=None)
data = pd.read_csv(os.path.join(args.output_prefix,"moods_filtered_scored_thresholded_clustered.bed"), sep="\t", header=None)
bw = pyBigWig.open(args.bigwig)
genome = pyfaidx.Fasta(args.genome)

print(data.shape)

corrs=[]
for i,r in tqdm(data.iterrows()):
	chr = r[0]
	start = r[1]
	end = r[2]
	key = "_".join(r[3].split("_")[0:2])

	val = np.nan_to_num(bw.values(chr,start,end)).reshape((-1,1))
	seq = one_hot.dna_to_one_hot([str(genome[chr][start:end])])[0]

	#score = np.abs(val)
	#trim_thresh = np.max(score) * 0.3
	#val[val<trim_thresh] = 0.0

	output = val*seq

	corr1 = np.max(signal.correlate2d(pfms[key][::-1, ::-1],output, mode="valid"))
	corr2 = np.max(signal.correlate2d(pfms[key],output, mode="valid"))

	corr = np.max([corr1, corr2])
	corrs.append(corr/(end-start))
	#corrs.append(corr)
	#print(corrs[-1])

#print(data.head())
#print(len(corrs))

data[9] = corrs
#print(data.shape)

#data.to_csv(os.path.join(args.output_prefix,"moods_filtered_scored_thresholded_clustered_cwm_actvations_normed_length.bed"),sep="\t", header=False, index=False)
#data.to_csv(os.path.join(args.output_prefix,"moods_filtered_scored_thresholded_clustered_cwm_actvations_subset_normed.bed"),sep="\t", header=False, index=False)
data.to_csv(os.path.join(args.output_prefix,"moods_filtered_scored_thresholded_clustered_cwm_actvations_normed.bed"),sep="\t", header=False, index=False)


new_data = []

clusters = list(set(data[8].values))


for idx in tqdm(clusters):
	obc = data[data[8]==idx].reset_index(drop=True)

	if len(obc) == 1:
		new_data.append(obc.loc[0,:].values)
	elif len(obc) == 0:
		continue
	else:
		while len(obc) > 0:
			idx = obc[9].idxmax()
			new_data.append(obc.loc[idx,:].values)
			width = (obc.loc[idx,2] - obc.loc[idx,1])//4
			obc = obc[ (obc[2] < (obc.loc[idx,1]+width)) | (obc[1] > (obc.loc[idx,2]-width))].reset_index(drop=True)

new_data = pd.DataFrame(new_data)

#new_data.to_csv(os.path.join(args.output_prefix,"overlaps_resolved_based_on_cwm_activations_normed_length_normed.bed"),sep="\t", header=False, index=False)
#new_data.to_csv(os.path.join(args.output_prefix,"overlaps_resolved_based_on_cwm_activations_subset_normed.bed"),sep="\t", header=False, index=False)
new_data.to_csv(os.path.join(args.output_prefix,"overlaps_resolved_based_on_cwm_activations_normed.bed"),sep="\t", header=False, index=False)





