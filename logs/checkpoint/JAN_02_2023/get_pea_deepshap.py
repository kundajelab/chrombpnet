
import pandas as pd
import os
import deepdish as dd
import numpy as np
from tqdm import tqdm
import pyfaidx
import one_hot

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

data = pd.read_csv("model_dir_atac.csv",header=None)
ddtpe="ATAC"
ddtpen=ddtpe+"_PE"
#cell_types=["HEPG2", "K562", "GM12878", "H1ESC", "IMR90"]
#cell_types=[ "H1ESC", "IMR90"]
#cell_types=["K562", "GM12878", "H1ESC", "IMR90", "HEPG2"]
#cell_types=["IMR90"]
#itype="counts"

#data = pd.read_csv("model_dir_dnase.csv",header=None)
#data = pd.read_csv("v1/model_dir_dnase_v2_interpret.csv",header=None)
#ddtpe="DNASE"
#ddtpen=ddtpe+"_SE"
#ddtpen=ddtpe+"_SE"
#cell_types=["HEPG2", "K562", "GM12878", "H1ESC", "IMR90"]
#cell_types=["K562", "GM12878", "H1ESC", "IMR90", "HEPG2"]
#cell_types=["HEPG2", "K562"]
#cell_types=["K562"]
#cell_types=["GM12878", "IMR90"]
#cell_types=["GM12878_new", "IMR90_new", "H1ESC_new"]
cell_types=["GM12878", "IMR90", "H1ESC"]
#cell_types=[ "IMR90"]
itype="counts"



#data = pd.read_csv("model_dir_dnase.csv",header=None)

genome_fa="/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa"



def get_seq(peaks_df, genome, width):
    """
    Same as get_cts, but fetches sequence from a given genome.
    """
    vals = []
    peaks_used = []
    for i, r in peaks_df.iterrows():
        sequence = str(genome[r[0]][(r[1]+r[9] - width//2):(r[1] + r[9] + width//2)])
        if len(sequence) == width:
            vals.append(sequence)
            peaks_used.append(True)
        else:
            peaks_used.append(False)
    assert(len(peaks_used)==peaks_df.shape[0])
    return one_hot.dna_to_one_hot(vals)


def filter_regions_to_peaks(bed_of_interest, merged, scores, output_prefix):

	output_prefix=output_prefix+'/'+cell_type+'_ATAC_in_peaks'
	print(output_prefix)
	boi = bed_of_interest[["chr", "start", "end", "summit"]].to_numpy().tolist()
	merged_val = merged[[0,1,2,9]].to_numpy().tolist()
	
	indices=[]
	dups = []
	#for i, val in enumerate(tqdm(merged_val)):
	for i, val in enumerate(merged_val):
		if val in boi:
			if val not in dups:
				indices.append(i)
				dups.append(val)

	print(len(indices))
	print(len(merged_val))
	print(len(boi))
	#assert(len(indices)==len(boi))
	merged.iloc[indices].to_csv(output_prefix+"."+itype+".interpreted_regions.bed", header=False, index=False, sep="\t")

	sub_scores = {
			'shap': {'seq': scores[indices]},
		}

	print(sub_scores['shap']['seq'].shape)

	dd.io.save(output_prefix+"."+itype+"_scores_new_compressed.h5",
					sub_scores,
					compression='blosc')
				
	
for cell_type in cell_types:
	ndata = data[data[1]==cell_type].reset_index()
	#cell_type = cell_type+"_new"
	bed_of_interest = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"+ddtpen+"/"+cell_type+"/data/peaks_no_blacklist.bed", sep="\t", header=None, names=NARROWPEAK_SCHEMA).astype(str)
	one_hots=None
	for i,r in ndata.iterrows():
		print(i,r[2])
		
		beds_path = os.path.join(r[2],"chrombpnet_model/interpret/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		if os.path.exists(beds_path):
			beds = pd.read_csv(beds_path, sep="\t", header=None)
		elif os.path.exists(os.path.join(r[2],"chrombpnet_model/interpret/merged."+cell_type+".interpreted_regions.bed")):
			beds_path = os.path.join(r[2],"chrombpnet_model/interpret/merged."+cell_type+".interpreted_regions.bed")
			beds = pd.read_csv(beds_path, sep="\t", header=None)
		else:
			beds_path = os.path.join(r[2],"interpret/merged."+cell_type+".interpreted_regions.bed")
			beds = pd.read_csv(beds_path, sep="\t", header=None)

		print(beds.head())
		beds["key"] = beds[0] + "_" + beds[1].astype(str) + "_" + beds[2].astype(str) + "_" + + beds[9].astype(str)

		ppath = os.path.join(r[2],"chrombpnet_model/interpret/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		if os.path.exists(ppath):
			scores = dd.io.load(ppath)
		elif os.path.exists(os.path.join(r[2],"interpret/merged."+cell_type+"."+itype+"_scores_new_compressed.h5")):
			ppath = os.path.join(r[2],"interpret/merged."+cell_type+"."+itype+"_scores_new_compressed.h5")
			scores = dd.io.load(ppath)
		else:
			ppath = os.path.join(r[2],"chrombpnet_model/interpret/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
			scores = dd.io.load(ppath)
	
		if i == 0 :
			output = scores['shap']['seq']
			shapez = output.shape
			init_beds = beds
			print(scores.keys())
			#raw = scores['raw']['seq']
			if 'raw' in scores:
				if one_hots is None:
					one_hots = scores['raw']['seq']
					print("one hots found")
		else:
			assert((init_beds["key"].to_numpy() == beds["key"].to_numpy()).all())
			print(scores['shap']['seq'].shape)
			assert(shapez==scores['shap']['seq'].shape)
			#assert(raw==scores['raw']['seq'])
			if 'raw' in scores:
				if one_hots is None:
					one_hots = scores['raw']['seq']
					print("one hots found")

			output = scores['shap']['seq']

		print(output.shape)

		os.makedirs(os.path.join(r[2],"interpret_peaks"), exist_ok=True)
		print(os.path.join(r[2],"interpret_peaks"))

		filter_regions_to_peaks(bed_of_interest, beds.astype(str), scores['shap']['seq'], os.path.join(r[2],"interpret_peaks"))

		del scores
