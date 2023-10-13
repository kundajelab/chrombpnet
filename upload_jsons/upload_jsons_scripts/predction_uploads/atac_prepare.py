import os
import json

encids = open("/mnt/lab_data2/anusri/chrombpnet/upload_jsons/upload_jsons_scripts/model_uploads/chrombpnet_models/data/atac_passed.txt").readlines()
encids = [line.strip() for line in encids]

odir="atac_temp/"


ignore_list = []
ignore_list += ["ENCSR955JSO"] # no .png
ignore_list += ["ENCSR459AYS"] # not all data paths present
ignore_list += ["ENCSR212LYK", "ENCSR260ZIV"] # fold log batch, train_chrombpnet_log
ignore_list += ["ENCSR637XSC", "ENCSR291GJU", "ENCSR200OML", "ENCSR868FGK"]


for encid in encids:
	ofile = "/scratch/users/anusri/chromatin-atlas-2022/DNASE/"+encid+"/prediction_uploads/average_preds/"+encid+".mean_preds_wo_bias.stat"
	if encid in ignore_list:
		continue
	if os.path.isfile(ofile):
		wbias = "/scratch/users/anusri/chromatin-atlas-2022/DNASE/"+encid+"/prediction_uploads/average_preds/"+encid+".mean_preds_w_bias.bw"
		nobias = "/scratch/users/anusri/chromatin-atlas-2022/DNASE/"+encid+"/prediction_uploads/average_preds/"+encid+".mean_preds_wo_bias.bw"

		assert(os.path.isfile(wbias)==True)	
		assert(os.path.isfile(nobias)==True)	

		output_json = {}
		output_json["experiment"] = encid
		output_json["predicted signal profile bigWig"] = wbias
		output_json["bias-corrected predicted signal profile bigWig"] = nobias


