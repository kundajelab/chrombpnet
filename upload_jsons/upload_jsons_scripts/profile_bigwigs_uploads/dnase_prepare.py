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

odir="temp/"
#odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/DNASE/production_uplaods/profile_interpret_bws_uploads/"

for encid in encids:
	bw_path = os.path.join("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/", encid+"/interpret_uploads/average_deepshap/"+encid+".mean_fold_profile.bw" )
	if encid in ignore_list:
		continue

	if os.path.isfile(bw_path):
		output_json = {}
		output_json["experiment"] = encid
		output_json["profile sequence contribution scores bigWig"] = bw_path
			
		if not os.path.isfile(odir+encid+".json"):
			f = open(odir+encid+".json", "w")
			json.dump(output_json, f, indent=4)
			f.close()
	else:
		print(encid)
