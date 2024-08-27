
import pandas as pd
import os
import deepdish as dd
import numpy as np
from tqdm import tqdm
import pyfaidx
import one_hot

#data = pd.read_csv("model_dir_atac.csv",header=None)
#ddtpe="ATAC"
#ddtpen=ddtpe+"_PE"
#cell_types=["HEPG2", "K562", "GM12878", "H1ESC", "IMR90"]
#cell_types=[ "H1ESC", "IMR90"]
#cell_types=["K562", "GM12878", "H1ESC", "IMR90", "HEPG2"]
#cell_types=["IMR90"]
#itype="counts"

#data = pd.read_csv("model_dir_dnase.csv",header=None)
data = pd.read_csv("v1/model_dir_dnase.csv",header=None)
ddtpe="DNASE"
ddtpen=ddtpe+"_PE"
#ddtpen=ddtpe+"_SE"
#cell_types=["HEPG2", "K562", "GM12878", "H1ESC", "IMR90"]
#cell_types=["K562", "GM12878", "H1ESC", "IMR90", "HEPG2"]
cell_types=["HEPG2", "K562"]
#cell_types=["K562"]
#cell_types=["GM12878", "IMR90"]
#cell_types=["GM12878_new", "IMR90_new", "H1ESC_new"]
#cell_types=["GM12878"]
#cell_types=["GM12878", "IMR90", "H1ESC"]
#cell_types=["H1ESC"]

itype="counts"



#data = pd.read_csv("model_dir_dnase.csv",header=None)

genome_fa="/mnt/lab_data2/anusri/chrombpnet/reference/hg38.genome.fa"

NARROWPEAK_SCHEMA=["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]

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
    #assert(le(peaks_used)==peaks_df.shape[0])
    print(sum(peaks_used), len(peaks_used))
    return one_hot.dna_to_one_hot(vals)


				
	
for cell_type in cell_types:
	ndata = data[data[1]==cell_type].reset_index()
	#cell_type = cell_type+"_new"
	one_hots=None
	for i,r in ndata.iterrows():
		print(i,r[2])

		
		beds_path = os.path.join(r[2],"chrombpnet_model/interpret/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		if os.path.exists(beds_path):
			try:
				beds = pd.read_csv(beds_path, sep="\t", header=None, compression='gzip')
			except:
				beds = pd.read_csv(beds_path, sep="\t", header=None)
				
		elif os.path.exists(os.path.join(r[2],"chrombpnet_model/interpret/merged."+cell_type+".interpreted_regions.bed")):
			beds_path = os.path.join(r[2],"chrombpnet_model/interpret/merged."+cell_type+".interpreted_regions.bed")
			beds = pd.read_csv(beds_path, sep="\t", header=None)
		else:
			beds_path = os.path.join(r[2],"interpret/merged."+cell_type+".interpreted_regions.bed")
			beds = pd.read_csv(beds_path, sep="\t", header=None)

		print(beds.shape)

		beds_path = os.path.join(r[2],"chrombpnet_model/interpret_ccre/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		if os.path.exists(beds_path):
			beds2 = pd.read_csv(beds_path, sep="\t", header=None)
		elif os.path.exists(os.path.join(r[2],"chrombpnet_model/interpret_ccre/merged."+cell_type+".interpreted_regions.bed")):
			beds_path = os.path.join(r[2],"chrombpnet_model/interpret_ccre/merged."+cell_type+".interpreted_regions.bed")
			beds2 = pd.read_csv(beds_path, sep="\t", header=None)
		else:
			irpath="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/"+r[2].split("/")[-1]
			beds_path = os.path.join(irpath,"chrombpnet_model/interpret_ccre/full_"+cell_type+".interpreted_regions_"+itype+".bed")
			beds2 = pd.read_csv(beds_path, sep="\t", header=None)

		print(beds2.shape)

		bedsf = pd.concat([beds, beds2])
		print(bedsf.shape)


		ppath = os.path.join(r[2],"chrombpnet_model/interpret/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		if os.path.exists(ppath):
			scores = dd.io.load(ppath)
		elif os.path.exists(os.path.join(r[2],"interpret/merged."+cell_type+"."+itype+"_scores_new_compressed.h5")):
			ppath = os.path.join(r[2],"interpret/merged."+cell_type+"."+itype+"_scores_new_compressed.h5")
			scores = dd.io.load(ppath)
		else:
			ppath = os.path.join(r[2],"interpret/merged."+cell_type+"."+itype+"_scores_new_compressed.h5")
			scores = dd.io.load(ppath)
	
		ppath = os.path.join(r[2],"chrombpnet_model/interpret_ccre/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		if os.path.exists(ppath):
			scores2 = dd.io.load(ppath)
		elif os.path.exists(os.path.join(r[2],"interpret_ccre/merged."+cell_type+"."+itype+"_scores_new_compressed.h5")):
			ppath = os.path.join(r[2],"interpret_ccre/merged."+cell_type+"."+itype+"_scores_new_compressed.h5")
			scores2 = dd.io.load(ppath)
		else:
			irpath="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/"+r[2].split("/")[-1]
			ppath = os.path.join(irpath,"chrombpnet_model/interpret_ccre/full_"+cell_type+"."+itype+"_scores.h5")			
			scores2 = dd.io.load(ppath)
	
		print(scores['shap']['seq'].shape)
		print(scores2['shap']['seq'].shape)


		scoresf = np.concatenate((scores['shap']['seq'], scores2['shap']['seq']), axis=0)
		print(scoresf.shape)

		assert(scoresf.shape[0]==bedsf.shape[0])

		counts_scores_dict = {
			'shap': {'seq': scoresf},
		}

		irpath="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_type+"/"+r[2].split("/")[-1]
		os.makedirs(os.path.join(irpath,"chrombpnet_model/interpret_all_with_ccre/"), exist_ok=True)
		print(os.path.join(irpath,"chrombpnet_model/interpret_all_with_ccre/"))
		output_prefix=os.path.join(irpath,"chrombpnet_model/interpret_all_with_ccre/full_"+cell_type+"."+itype+"_scores_new_compressed.h5")
		dd.io.save(output_prefix,
			counts_scores_dict,
			compression='blosc')

		output_prefix_bed=os.path.join(irpath,"chrombpnet_model/interpret_all_with_ccre/full_"+cell_type+".interpreted_regions_"+itype+".bed")
		bedsf.to_csv(output_prefix_bed, sep="\t", header=False, index=False, compression='gzip')
		#break
