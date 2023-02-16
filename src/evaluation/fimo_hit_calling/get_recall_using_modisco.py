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
            print(len(metacluster["seqlets"].value))
            if "patterns" not in metacluster["seqlets_to_patterns_result"]:
                continue
            patterns = metacluster["seqlets_to_patterns_result"]["patterns"]
            num_patterns = len(patterns["all_pattern_names"][:])
            for pattern_i, pattern_name in enumerate(patterns["all_pattern_names"][:]):
                pattern_name = pattern_name
                pattern = patterns[pattern_name]
                key = "metacluster_"+str(metacluster_i)+".pattern_"+str(pattern_i)
                
                if key not in pfms:
                    pfms[key] = []
                    
                for seqlet in pattern["seqlets_and_alnmts"]["seqlets"]:
                    pfms[key].append(seqlet)
               
    return pfms
    
    

def main():
	args = fetch_arguments()      

	tfm_results_path=args.modisco_obj
	pfms = import_tfmodisco_motifs(tfm_results_path)
	tomtom=args.tomtom_anotation_path
	tomtom = pd.read_csv(tomtom, sep="\t", header=0)
	bed = args.regions
	bed = pd.read_csv(bed, sep="\t", header=None)

	## map seqlets to actual coordinates 
	lists=[]
	for key in pfms:
		print(key)
		if key in tomtom["Pattern"].values:
			match_name = tomtom[tomtom["Pattern"]==key]["Match_1"]
			keyd = key.split("_")[1].replace(".pattern","")+"_"+key.split("_")[-1]+"_"+match_name
			keyd = keyd.values[0]
			#print(keyd.values[0])
			for seqlet in pfms[key]:
				vals = seqlet.split(",")
				peak_id = int(vals[0].split(":")[1])
				ss = vals[1].split(":")[1]
				ee = vals[2].split(":")[1]
				blist = [bed.loc[peak_id,0], bed.loc[peak_id,1]+bed.loc[peak_id,9]-250+int(ss), bed.loc[peak_id,1]+bed.loc[peak_id,9]-250+int(ee), keyd]
				lists.append(blist)
		else:
			print(key)
        
	labels = pd.DataFrame(lists)
	labels.to_csv(os.path.join(args.output_dir,"modisco_annotations.bed"), sep="\t", header=False, index=False)
	labels.columns = ["tfchr", "tfstart", "tfend", "tfname"]

	f = open(os.path.join(args.output_dir,"hit_calls_in_modisco_annotations.bed"), "w")
	comm = ["bedtools", "intersect", "-f", "1.0", "-a"]
	comm += [args.hit_calls_bed]
	comm += ["-b"]
	comm += [os.path.join(args.output_dir,"modisco_annotations.bed")]
	comm += ["-wo"]
	with open(os.path.join(args.output_dir,"hit_calls_in_modisco_annotations.bed"), "w") as f:
		proc = subprocess.Popen(comm, stdout=f)
		proc.wait()
	
	
	names = ["hit_chr", "hit_start", "hit_end", "hit_name", "strand", "imp_score", "cwm_score", "tfchr", "tfstart", "tfend", "tfname", "nb"]
	hits = pd.read_csv(os.path.join(args.output_dir,"hit_calls_in_modisco_annotations.bed"), sep="\t", header=None, names=names).drop_duplicates()
	all_uniq = np.array(list(set(hits["hit_name"].values.tolist())))
	pattern_ids = [int(x.split("_")[1]) for x in all_uniq]
	indxs = np.argsort(pattern_ids)
	


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
	
		dwms = temp[temp["hit_name"]==motif].drop_duplicates(subset=["tfchr", "tfstart", "tfend", "tfname"])
		dwms["id"] = dwms["tfchr"].astype(str)+"_"+dwms["tfstart"].astype(str)+"_"+dwms["tfend"].astype(str)


		recovered_hits = dwms.shape[0]
   
		# all failed hits 
		print("failure cases example for motif: "+motif)
	
		ttms = temp[~temp["id"].isin(dwms["id"])].drop_duplicates(subset=["tfchr", "tfstart", "tfend", "tfname"])
		print(ttms.head(100))
	
		f.write(",".join([motif,str(np.round(recovered_hits*100/modisco_hits.shape[0],2))]))
		f.write("\n")

if __name__=="__main__":
	main()