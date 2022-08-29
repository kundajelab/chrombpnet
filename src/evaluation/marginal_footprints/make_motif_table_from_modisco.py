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
modisco_temp = '/oak/stanford/groups/akundaje/projects/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/GM12878/modisco_crop_500_100K_seqs_1/'
footprint_url_text="http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/GM12878/footprints/"

celltypes=["HEPG2", "K562", "IMR90", "H1ESC", "GM12878"]

default_chars = {0:"A", 1:"C", 2:"G", 3:"T"}

rows = []

for idx,celltype_main in enumerate(celltypes):

	celltype = celltypes[idx]+"_COUNTS"
	url_text = "http://mitra.stanford.edu/kundaje/anusri/chrombpnet_paper_new/modisco_jun_30/modisco/ATAC/"+celltypes[idx]+"/modisco_crop_500_100K_seqs_1"
	modisco_dir=modisco_temp.replace("GM12878",celltypes[idx])

	tomtom_path=os.path.join(modisco_dir,"counts.tomtom.tsv")
	tomtom=pd.read_csv(tomtom_path, sep="\t")
	label_dict = {}
	for index,row in tomtom.iterrows():
		label_dict[row['Pattern']] = row['Match_1']
	modisco_results = h5py.File(os.path.join(modisco_dir,'modisco_results_allChroms_counts.hdf5'), 'r')

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
			logo1="=IMAGE("+"\""+url_text+"/modisco_logos/counts."+pattern_name+".cwm.fwd.png"+"\""+",1)"
			logo2="=IMAGE("+"\""+url_text+"/modisco_logos/counts."+pattern_name+".cwm.rev.png"+"\""+",1)"

			gm_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "GM12878")+"/GM12878."+row_name+".footprint.png"+"\""+",1)"
			hep_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "HEPG2")+"/HEPG2."+row_name+".footprint.png"+"\""+",1)"
			k562_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "K562")+"/K562."+row_name+".footprint.png"+"\""+",1)"
			h1_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "H1ESC")+"/H1ESC."+row_name+".footprint.png"+"\""+",1)"
			imr90_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "IMR90")+"/IMR90."+row_name+".footprint.png"+"\""+",1)"

			dnase_gm_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "GM12878").replace("ATAC","DNASE")+"/GM12878."+row_name+".footprint.png"+"\""+",1)"
			dnase_hep_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "HEPG2").replace("ATAC","DNASE")+"/HEPG2."+row_name+".footprint.png"+"\""+",1)"
			dnase_k562_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "K562").replace("ATAC","DNASE")+"/K562."+row_name+".footprint.png"+"\""+",1)"
			dnase_h1_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "H1ESC").replace("ATAC","DNASE")+"/H1ESC."+row_name+".footprint.png"+"\""+",1)"
			dnase_imr90_footprint="=IMAGE("+"\""+footprint_url_text.replace("GM12878", "IMR90").replace("ATAC","DNASE")+"/IMR90."+row_name+".footprint.png"+"\""+",1)"


			rows.append([row_name, "".join(seqs), str(label_dict[metacluster_name+"."+pattern_name]), logo1, logo2, gm_footprint, 0 ,dnase_gm_footprint, 0,hep_footprint,0,dnase_hep_footprint,0, k562_footprint,0,dnase_k562_footprint,0, imr90_footprint,0,dnase_imr90_footprint,0, h1_footprint,0,dnase_h1_footprint,0])
			print(rows[-1])


df = pd.DataFrame(rows, columns=["MOTIF_ID", "MOTIF_SEQ","BEST_TOMTOM_MATCH" ,"CWM_FWD", "CWM_REV", "GM12878_ATAC_FOOTPRINT","GM12878_ATAC_FOOTPRINT_SCORE", "GM12878_DNASE_FOOTPRINT","GM12878_DNASE_FOOTPRINT_SCORE","HEPG2_ATAC_FOOTPRINT","HEPG2_ATAC_FOOTPRINT_SCORE","HEPG2_DNASE_FOOTPRINT","HEPG2_DNASE_FOOTPRINT_SCORE", "K562_ATAC_FOOTPRINT","K562_ATAC_FOOTPRINT_SCORE","K562_DNASE_FOOTPRINT","K562_DNASE_FOOTPRINT_SCORE","IMR90_ATAC_FOOTPRINT","IMR90_ATAC_FOOTPRINT_SCORE","IMR90_DNASE_FOOTPRINT","IMR90_DNASE_FOOTPRINT_SCORE", "H1_ATAC_FOOTPRINT","H1_ATAC_FOOTPRINT_SCORE", "H1_DNASE_FOOTPRINT","H1_DNASE_FOOTPRINT_SCORE"])
df.to_csv("union_motifs_1.tsv", sep="\t", header=False, index=False)

