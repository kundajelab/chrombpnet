import os
import json

primary_encids = open("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_jsons_scripts/model_uploads/chrombpnet_models/data/primary_passed.txt").readlines()
primary_encids = [line.strip() for line in primary_encids]

celline_encids = open("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_jsons_scripts/model_uploads/chrombpnet_models/data/cellline_passed.txt").readlines()
celline_encids = [line.strip() for line in celline_encids]

tissue_encids = open("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_jsons_scripts/model_uploads/chrombpnet_models/data/tissue_passed.txt").readlines()
tissue_encids = [line.strip() for line in tissue_encids]

invitro_encids = open("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_jsons_scripts/model_uploads/chrombpnet_models/data/invitro_passed.txt").readlines()
invitro_encids = [line.strip() for line in invitro_encids]

encids = tissue_encids + primary_encids + celline_encids + invitro_encids
#encids = primary_encids + celline_encids + tissue_encids

ignore_list = []
# missing args file
ignore_list += ["ENCSR720TCN", "ENCSR241BNZ", "ENCSR516JCM", "ENCSR642DZF", "ENCSR787ERP"]
# missing train_test
ignore_list += ['ENCSR000EOT']
# chrombpnet models
ignore_list += ["ENCSR000EMT", "ENCSR149XIL", "ENCSR477RTP", "ENCSR000EOT", "ENCSR000EMU"]

#odir="temp/"
odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/DNASE/production_uplaods/counts_interpret_bw_uploads/sep_25_2023/"

def make_bb_file(in_bed, out_bb):
	assert(os.path.isfile("temp.bed")==False)
	command = "zcat "+in_bed+" | LC_COLLATE=C sort -k1,1 -k2,2n > temp.bed"
	print(command)
	os.system(command)

	command = "bedToBigBed temp.bed /oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/reference/chrom.sizes "+out_bb
	print(command)
	os.system(command)


	command = "rm temp.bed"
	print(command)
	os.system(command)

for encid in encids:
	bw_path = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/", encid+"/interpret_uploads/average_deepshap/"+encid+".mean_fold_counts.bw" )
	sel_path = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/", encid+"/interpret_uploads/selected_new.regions.valid.merged.bed.gz" )
	sel_path_bb = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/", encid+"/interpret_uploads/selected_new.regions.valid.merged.bigBed" )
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
		print(encid)
