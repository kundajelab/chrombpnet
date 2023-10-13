import h5py 
import numpy as np
import pandas as pd
import os
import argparse
import subprocess

def fetch_arguments():
	parser=argparse.ArgumentParser(description="get recall with modisco")
	parser.add_argument("-mo", "--modisco-obj", type=str, required=True, help="Modisco object to extract pfms")
	parser.add_argument("-t", "--tomtom_anotation_path", type=str, required=False, default=None, help="If provided hits are renamed with tomtom annotations instead of modisco annotations")
	parser.add_argument("-r", "--regions", type=str, required=True, help="10 column bed file of peaks. The order of the peaks should correspond to the order of contribution scores regions provided for TFModisco.")
	parser.add_argument("-ht", "--hit_calls_bed", type=str, required=True, help="Output of hitcalling.py")
	parser.add_argument("-o", "--output-dir", type=str, required=True, help="Output directory to store results")
	args = parser.parse_args()
	return args
	
def import_tfmodisco_motifs(tfm_results_path, trim=True, only_pos=True):
    pfms = {}
    with h5py.File(tfm_results_path, "r") as f:
        metaclusters = f["metacluster_idx_to_submetacluster_results"]
        num_metaclusters = len(metaclusters.keys())
        for metacluster_i, metacluster_key in enumerate(metaclusters.keys()):
            metacluster = metaclusters[metacluster_key]
            #print(len(metacluster["seqlets"].value))
            if "patterns" not in metacluster["seqlets_to_patterns_result"]:
                continue
            patterns = metacluster["seqlets_to_patterns_result"]["patterns"]
            num_patterns = len(patterns["all_pattern_names"][:])
            for pattern_i, pattern_name in enumerate(patterns["all_pattern_names"][:]):
                pattern_name = pattern_name
                pattern = patterns[pattern_name]
                key = "metacluster_"+str(metacluster_i)+".pattern_"+str(pattern_i)
                
                #if metacluster_i!=0 or pattern_i!=0:
                #b    continue
					               
                if key not in pfms:
                    pfms[key] = []
                    
                for seqlet in pattern["seqlets_and_alnmts"]["seqlets"]:
                    #pfms[key].append(seqlet.decode("utf-8"))
                    pfms[key].append(seqlet)
                    
                #break
               
    return pfms
    
    

def main():
	args = fetch_arguments()      

	tfm_results_path=args.modisco_obj
	pfms = import_tfmodisco_motifs(tfm_results_path)
	tomtomf=args.tomtom_anotation_path
	if tomtomf is not None:
		tomtom = pd.read_csv(tomtomf, sep="\t", header=0)
	else: 
		tomtom = None
	bed = args.regions
	bed = pd.read_csv(bed, sep="\t", header=None)

	## map seqlets to actual coordinates 
	lists=[]
	for key in pfms:
		print(key)
		if tomtom is not None:
			#key in tomtom["Pattern"].values
			match_name = tomtom[tomtom["Pattern"]==key]["Match_1"]
			keyd = key.split("_")[1].replace(".pattern","")+"_"+key.split("_")[-1]+"_"+match_name
			print(keyd)
			keyd = keyd.values[0]
			#print(keyd.values[0])
			for seqlet in pfms[key]:
				vals = seqlet.split(",")
				peak_id = int(vals[0].split(":")[1])
				ss = vals[1].split(":")[1]
				ee = vals[2].split(":")[1]
				blist = [bed.loc[peak_id,0], bed.loc[peak_id,1]+bed.loc[peak_id,9]-500+int(ss), bed.loc[peak_id,1]+bed.loc[peak_id,9]-500+int(ee), keyd]
				lists.append(blist)
		else:
			keyd = key.split("_")[1].replace(".pattern","")+"_"+key.split("_")[-1]
			#keyd = keyd.values[0]
			#print(keyd)
			for seqlet in pfms[key]:
				#print(seqlet)
				vals = seqlet.split(",")
				peak_id = int(vals[0].split(":")[1])
				ss = vals[1].split(":")[1]
				ee = vals[2].split(":")[1]
				blist = [bed.loc[peak_id,0], bed.loc[peak_id,1]+bed.loc[peak_id,9]-500+int(ss), bed.loc[peak_id,1]+bed.loc[peak_id,9]-500+int(ee), keyd]
				lists.append(blist)
			print(key)
        
	labels = pd.DataFrame(lists)
	labels.to_csv(os.path.join(args.output_dir,"modisco_annotations.bed"), sep="\t", header=False, index=False)
	labels.columns = ["tfchr", "tfstart", "tfend", "tfname"]

	f = open(os.path.join(args.output_dir,"hit_calls_in_modisco_annotations.bed"), "w")
	comm = ["bedtools", "intersect", "-f", "0.5", "-a"]
	comm += [args.hit_calls_bed]
	comm += ["-b"]
	comm += [os.path.join(args.output_dir,"modisco_annotations.bed")]
	comm += ["-wo"]
	with open(os.path.join(args.output_dir,"hit_calls_in_modisco_annotations.bed"), "w") as f:
		proc = subprocess.Popen(comm, stdout=f)
		proc.wait()

#	names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "match_score", "imp_score", "cwm_score", "pvalue", "qvalue", "cluster", "tfchr", "tfstart", "tfend", "tfname", "nb"]
#	names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "match_score", "imp_score", "cwm_score", "pvalue", "qvalue", "cluster" , "tfchr", "tfstart", "tfend", "tfname", "nb"]
	names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "match_score", "imp_score", "cwm_score", "pvalue", "qvalue", "tfchr", "tfstart", "tfend", "tfname", "nb"]
	
	#if tomtom is not None:
	#	names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "match_score", "imp_score", "cwm_score", "pvalue", "qvalue", "tfchr", "tfstart", "tfend", "tfname", "nb"]

	#else:
	#	#names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "imp_score", "peak_index", "cwm_score", "pvalue",  "qvalue","tfchr", "tfstart", "tfend", "tfname", "nb"]
	#	#names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "imp_score", "tfchr", "tfstart", "tfend", "tfname", "nb"]
	#	names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "imp_score", "peak_index", "cwm_score", "tfchr", "tfstart", "tfend", "tfname", "nb"]
	#	#names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "imp_score", "cwm_score", "tfchr", "tfstart", "tfend", "tfname", "nb"]

	hits = pd.read_csv(os.path.join(args.output_dir,"hit_calls_in_modisco_annotations.bed"), sep="\t", header=None, names=names).drop_duplicates()
	all_uniq = np.array(list(set(hits["hit_name"].values.tolist())))
	print(all_uniq)
	print(hits.head())
	pattern_ids = [int(x.split("_")[1]) for x in all_uniq]
	indxs = np.argsort(pattern_ids)
	
	motifs_names = []
	for motif in all_uniq[indxs]:
		
		if int(motif.split("_")[0]) != 0:
			continue
		motifs_names.append(motif)
	
	confusion_matrix = np.zeros((len(motifs_names), len(motifs_names)))


	motifs = []
	pwms_all=[]
	cwms_all=[]
	f = open(os.path.join(args.output_dir,"recall_per_motif.csv"),"w")
	for motif in all_uniq[indxs]:
		
		if int(motif.split("_")[0]) != 0:
			continue
		
		modisco_hits = labels[labels["tfname"]==motif].drop_duplicates()
	
		temp = pd.merge(modisco_hits, hits, left_on=["tfchr", "tfstart", "tfend", "tfname"],right_on=["tfchr", "tfstart", "tfend", "tfname"], how="inner").drop_duplicates()

		temp["id"] = temp["tfchr"].astype(str)+"_"+temp["tfstart"].astype(str)+"_"+temp["tfend"].astype(str)
		
		if motif == "0_3_BACH2_HUMAN.H11MO.0.A":
			temp["hit_name"].replace('0_28_ATF3_MOUSE.H11MO.0.A','0_3_BACH2_HUMAN.H11MO.0.A', inplace=True)

		if motif == "0_4_NFIA_HUMAN.H11MO.0.C":
			temp["hit_name"].replace('0_35_NFIC_HUMAN.H11MO.0.A','0_4_NFIA_HUMAN.H11MO.0.C', inplace=True)
		
		if motif == "0_35_NFIC_HUMAN.H11MO.0.A":
			temp["hit_name"].replace('0_4_NFIA_HUMAN.H11MO.0.C', '0_35_NFIC_HUMAN.H11MO.0.A', inplace=True)
		
		if motif == "0_5_Gabpa_MA0062.2":
			temp["hit_name"].replace('0_11_ETV4_MOUSE.H11MO.0.B', '0_5_Gabpa_MA0062.2', inplace=True)

		if motif == "0_11_ETV4_MOUSE.H11MO.0.B":
			temp["hit_name"].replace( '0_5_Gabpa_MA0062.2', '0_11_ETV4_MOUSE.H11MO.0.B', inplace=True)
		
		if motif == "0_9_ZN143_MOUSE.H11MO.0.A":
			temp["hit_name"].replace('0_17_ZNF76_HUMAN.H11MO.0.C', '0_9_ZN143_MOUSE.H11MO.0.A', inplace=True)
		
		if motif == "0_17_ZNF76_HUMAN.H11MO.0.C":
			temp["hit_name"].replace('0_9_ZN143_MOUSE.H11MO.0.A', '0_17_ZNF76_HUMAN.H11MO.0.C', inplace=True)
		
		if motif == "0_7_NRF1_HUMAN.H11MO.0.A":
			temp["hit_name"].replace('0_20_NRF1_MA0506.1', '0_7_NRF1_HUMAN.H11MO.0.A', inplace=True)
				
		if motif == "0_20_NRF1_MA0506":
			temp["hit_name"].replace('0_7_NRF1_HUMAN.H11MO.0.A', '0_20_NRF1_MA0506.1', inplace=True)
	
		if motif == "0_10_BCL6_HUMAN.H11MO.0.A":
			temp["hit_name"].replace('0_33_Stat5a+Stat5b_MA0519.1','0_10_BCL6_HUMAN.H11MO.0.A', inplace=True)	
		
		if motif == "0_33_Stat5a+Stat5b_MA0519.1":
			temp["hit_name"].replace('0_10_BCL6_HUMAN.H11MO.0.A', '0_33_Stat5a+Stat5b_MA0519.1', inplace=True)	
		
		dwms = temp[temp["hit_name"]==motif].drop_duplicates(subset=["tfchr", "tfstart", "tfend", "tfname"])
		dwms["id"] = dwms["tfchr"].astype(str)+"_"+dwms["tfstart"].astype(str)+"_"+dwms["tfend"].astype(str)


		recovered_hits = dwms.shape[0]
		print(recovered_hits)
   
		# all failed hits 
		print("failure cases example for motif: "+motif)
	
		ttms = temp[~temp["id"].isin(dwms["id"])].drop_duplicates(subset=["tfchr", "tfstart", "tfend", "tfname"])
		print(ttms.head(100)[["tfchr", "tfstart", "tfend", "tfname", "hit_chr", "hit_start", "hit_end", "hit_name"]])
		print(ttms["hit_name"].value_counts(dropna=False))
		
		if  motif == "0_2_KLF12_HUMAN.H11MO.0.C":
			ttms.to_csv(os.path.join(args.output_dir,"klf_versus_ctcf.bed"), sep="\t", header=True, index=False)
			
		if  motif == "0_2_KLF12_HUMAN.H11MO.0.C":
			ttms.to_csv(os.path.join(args.output_dir,"klf_versus_ctcf.bed"), sep="\t", header=True, index=False)
				
		if  motif == "0_2_KLF12_HUMAN.H11MO.0.C":
			dwms.to_csv(os.path.join(args.output_dir,"sp1.bed"), sep="\t", header=True, index=False)
		if  motif == "0_1_CTCF_MA0139.1":
			dwms.to_csv(os.path.join(args.output_dir,"ctcf.bed"), sep="\t", header=True, index=False)
		
		dicts = ttms["hit_name"].value_counts(dropna=False).to_dict()
		
		idx = motifs_names.index(motif)
		for key in dicts:
			if key in motifs_names:
				jdx = motifs_names.index(key)
			else:
				continue
			confusion_matrix[idx,jdx] = dicts[key]/modisco_hits.shape[0]
	
		f.write(",".join([motif,str(np.round(recovered_hits*100/modisco_hits.shape[0],2))]))
		f.write("\n")

	f = open(os.path.join(args.output_dir,"confusion_matrix.npy"), "wb")
	np.save(f,confusion_matrix)
	f = open(os.path.join(args.output_dir,"motif_names.npy"), "wb")
	np.save(f,motifs_names)
if __name__=="__main__":
	main()
