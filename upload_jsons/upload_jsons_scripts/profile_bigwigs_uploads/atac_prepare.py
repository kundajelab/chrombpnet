import os
import json

encids = open("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_jsons_scripts/model_uploads/chrombpnet_models/data/atac_passed.txt").readlines()
encids = [line.strip() for line in encids]

#odir="atac_temp/"
odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/ATAC/production_uploads/counts_interpret_bw_uploads/sep_25_2023/"

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

ignore_list = []
ignore_list += ["ENCSR955JSO"] # no .png
ignore_list += ["ENCSR459AYS"] # not all data paths present
ignore_list += ["ENCSR212LYK", "ENCSR260ZIV"] # fold log batch, train_chrombpnet_log
ignore_list += ["ENCSR637XSC", "ENCSR291GJU", "ENCSR200OML", "ENCSR868FGK"]


for encid in encids:
	if not os.path.isfile("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/ATAC/production_uploads/model_uploads/aug_28_2023/chrombpnet_model/"+encid+".json"):
		continue

	bw_path = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/", encid+"/interpret_uploads/average_deepshap/"+encid+".mean_fold_counts.bw" )
	sel_path = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/", encid+"/interpret_uploads/selected.regions.valid.merged.bed.gz" )
	sel_path_bb = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/", encid+"/interpret_uploads/selected.regions.valid.merged.bigBed" )

	if encid in ignore_list:
		continue

	if os.path.isfile(bw_path):
		output_json = {}
		output_json["experiment"] = encid
		output_json["counts sequence contribution scores bigWig"] = bw_path
		if os.path.isfile(sel_path):
			output_json["selected regions for predicted signal and sequence contribution scores bed"] = sel_path

			if os.path.isfile(sel_path_bb):
				output_json["selected regions for predicted signal and sequence contribution scores bigBed"] = sel_path_bb
			else:
				make_bb_file(sel_path, sel_path_bb)
				output_json["selected regions for predicted signal and sequence contribution scores bigBed"] = sel_path_bb
			
			if not os.path.isfile(odir+encid+".json"):
				f = open(odir+encid+".json", "w")
				json.dump(output_json, f, indent=4)
				f.close()

	else:
		if os.path.isfile("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/ATAC/production_uploads/model_uploads/aug_28_2023/chrombpnet_model/"+encid+".json"):
			print(encid)
