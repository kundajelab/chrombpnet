import numpy as np
import pandas as pd
import h5py
import pysam
import os
from modisco.visualization import viz_sequence
from modisco import util
from matplotlib import pyplot as plt
import pybedtools

pd.options.display.max_rows = 500
pd.options.display.max_columns = 500

mode="counts"
sourcetypes=["BIAS", "SIGNAL"]
rows = []

default_chars = {0:"A", 1:"C", 2:"G", 3:"T"}

for idx,srct in enumerate(sourcetypes):

	#tomtom_path=os.path.join(modisco_dir,"counts.tomtom.tsv")
	#tomtom=pd.read_csv(tomtom_path, sep="\t")
	#label_dict = {}
	#for index,row in tomtom.iterrows():
	#	label_dict[row['Pattern']] = row['Match_1']

	if srct == "SIGNAL":
		modisco_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/GM12878/GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0/'
		modisco_results = h5py.File(os.path.join(os.path.join(modisco_dir, srct),'modisco_crop_500/modisco_results_allChroms_counts.hdf5'), 'r')
		celltype = "GM12878_"+srct

	if srct == "BIAS":
		modisco_dir = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/ATAC_PE/K562/K562_02.17.2022_bias_128_4_1234_0.5_fold_0'
		modisco_results = h5py.File(os.path.join(os.path.join(modisco_dir, srct),'modisco/modisco_results_allChroms_profile.hdf5'), 'r')
		celltype = "K562_"+srct


	for metacluster_name in modisco_results["metacluster_idx_to_submetacluster_results"]:
		if metacluster_name == "metacluster_1":
			continue
		metacluster = modisco_results["metacluster_idx_to_submetacluster_results"][metacluster_name]
		all_pattern_names = [x.decode("utf-8") for x in list(metacluster["seqlets_to_patterns_result"]["patterns"]["all_pattern_names"][:])]

		for pattern_name in all_pattern_names:
			print(pattern_name)
			row_name=celltype+"_"+metacluster_name+"_"+pattern_name
			ppm = np.array(metacluster['seqlets_to_patterns_result']['patterns'][pattern_name]['sequence']['fwd'])
			cwm_fwd = np.array(metacluster['seqlets_to_patterns_result']['patterns'][pattern_name]["task0_contrib_scores"]['fwd'])
			score_fwd = np.sum(np.abs(cwm_fwd), axis=1)
			trim_thresh_fwd = np.max(score_fwd) * 0.3
			pass_inds_fwd = np.where(score_fwd >= trim_thresh_fwd)[0]
			start_fwd, end_fwd = max(np.min(pass_inds_fwd) - 4, 0), min(np.max(pass_inds_fwd) + 4 + 1, len(score_fwd) + 1)
			trimmed_cwm_fwd = cwm_fwd[start_fwd:end_fwd]

			seqs = np.vectorize(default_chars.get)(np.argmax(trimmed_cwm_fwd,axis=1))

			rows.append([row_name, "".join(seqs)])
			print(rows[-1])


df = pd.DataFrame(rows, columns=["MOTIF_ID", "MOTIF_SEQ"])
df.to_csv("gm_benchmarking_motifs.tsv", sep="\t", header=True, index=False)

