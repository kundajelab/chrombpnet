import os
import json


encids = ["K562", "HEPG2", "IMR90_new", "H1ESC_new", "GM12878_new"]

encode_id = {"HEPG2": "ENCSR149XIL",
        "K562": "ENCSR000EOT",
        "IMR90_new": "ENCSR477RTP",
        "GM12878_new": "ENCSR000EMT",
        "H1ESC_new": "ENCSR000EMU"}
odir='dnase/'
for encid in encids:
	ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/preds_upload/average_preds_with_ccre_vf/"+encid+".mean_preds_wo_bias.stat"
	if os.path.isfile(ofile):
		print(encid)
		wbias = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/preds_upload/average_preds_with_ccre_vf/"+encid+".mean_preds_w_bias.bw"
		nobias = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+encid+"/preds_upload/average_preds_with_ccre_vf/"+encid+".mean_preds_wo_bias.bw"

		assert(os.path.isfile(wbias)==True)	
		assert(os.path.isfile(nobias)==True)	

		output_json = {}
		output_json["experiment"] = encode_id[encid]
		output_json["predicted signal profile bigWig"] = wbias
		output_json["bias-corrected predicted signal profile bigWig"] = nobias
		
		if not os.path.isfile(odir+encode_id[encid]+".json"):
			f = open(odir+encode_id[encid]+".json", "w")
			json.dump(output_json, f, indent=4)
			f.close()
