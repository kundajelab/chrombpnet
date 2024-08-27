import os
import json
import pandas as pd
import pybedtools

#encids = ["K562", "HEPG2", "IMR90", "H1ESC", "GM12878"]
#encids = ["K562", "HEPG2"]
encids = ["K562", "HEPG2", "IMR90_new", "H1ESC_new", "GM12878_new"]

encode_id = {"HEPG2": "ENCSR149XIL",
        "K562": "ENCSR000EOT",
        "IMR90_new": "ENCSR477RTP",
        "GM12878_new": "ENCSR000EMT",
        "H1ESC_new": "ENCSR000EMU"}
odir='dnase/'
def make_bb_file(in_bed, out_bb):
	assert(os.path.isfile("atac_temp.bed")==False)
	command = "zcat "+in_bed+" | LC_COLLATE=C sort -k1,1 -k2,2n > atac_temp.bed"
	print(command)
	os.system(command)

	command = "bedToBigBed atac_temp.bed /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/reference/chrom.sizes "+out_bb
	print(command)
	os.system(command)


	command = "rm atac_temp.bed"
	print(command)
	os.system(command)

chrs = list(map(str,list(range(1,23))))
chrs = ['chr'+i for i in chrs]
chrs = chrs + ['chrX', 'chrY']

def make_sel_bedfile(in_bed, out_bed):
	input_bed = pd.read_csv(in_bed, sep="\t", header=None)
	print(input_bed.shape)
	input_bed = input_bed[input_bed[0].isin(chrs)]
	print(input_bed.shape)

	input_bed[1] = input_bed[1]+input_bed[9]-500
	input_bed[2] = input_bed[1]  + 1000
	print(input_bed.head())
	x = pybedtools.BedTool.from_dataframe(input_bed[[0,1,2]])
	x = x.sort().merge()
	output_bed = x.to_dataframe()
	print(output_bed.shape)
	print(output_bed.head())
	output_bed.to_csv(out_bed, sep='\t', header=False, index=False)

for encid in encids:
	print(encid)
	ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.counts_scores_new_compressed.stats"
	if os.path.isfile(ofile):
		counts_bw = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.counts_scores_new_compressed.bw"
	else:
		print(ofile)
		ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.counts_scores.stats"
		if os.path.isfile(ofile):
			counts_bw = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.counts_scores.bw"
		else:
			counts_bw = None
			print(ofile)
			continue

		
	ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores_new_compressed.stats"
	if os.path.isfile(ofile):
		profile_bw = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores_new_compressed.bw"
	else:	
		print(ofile)
		ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores.stats"
		if os.path.isfile(ofile):
			profile_bw = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores.bw"
		else:
			profile_bw = None
			print(ofile)

			continue
	
	ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/interpret_upload/average_preds/per_folds.inputs.bed.gz"
	if os.path.isfile(ofile):
		sel_path = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/", encid+"/interpret_upload/average_preds/selected.regions.valid.merged.bed.gz" )
		sel_path_bb = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/", encid+"/interpret_upload/average_preds/selected.regions.valid.merged.bigBed" )
		if not os.path.isfile(sel_path):
			make_sel_bedfile(ofile, sel_path)
		
		if os.path.isfile(sel_path) and (not os.path.isfile(sel_path_bb)):
			make_bb_file(sel_path, sel_path_bb)
		
	else:
		sel_path=None
		sel_path_bb=None
		print(ofile)
		continue

	assert(os.path.isfile(counts_bw)==True)	
	assert(os.path.isfile(profile_bw)==True)	
	

	output_json = {}
	output_json["experiment"] = encode_id[encid]
	output_json["counts sequence contribution scores bigWig"] = counts_bw
	output_json["profile sequence contribution scores bigWig"] = profile_bw
		
	if os.path.isfile(sel_path_bb):
		output_json["selected regions for predicted signal and sequence contribution scores bigBed"] = sel_path_bb

	if os.path.isfile(sel_path):
		output_json["selected regions for predicted signal and sequence contribution scores bed"] = sel_path

	if not os.path.isfile(odir+encode_id[encid]+".json"):
		f = open(odir+encode_id[encid]+".json", "w")
		json.dump(output_json, f, indent=4)
		f.close()
