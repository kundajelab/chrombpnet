import os
import json


encids = ["K562", "HEPG2", "IMR90", "H1ESC", "GM12878"]

encode_id = {"K562": "ENCSR868FGK",
"GM12878": "ENCSR637XSC",
"HEPG2": "ENCSR291GJU",
"IMR90": "ENCSR200OML",
"H1ESC": "GSE267154"}  

odir='atac/'
for encid in encids:
	ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+encid+"/preds_upload/average_preds/"+encid+".mean_preds_wo_bias.stat"
	if os.path.isfile(ofile):
		print(encid)
		wbias = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+encid+"/preds_upload/average_preds/"+encid+".mean_preds_w_bias.bw"
		nobias = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+encid+"/preds_upload/average_preds/"+encid+".mean_preds_wo_bias.bw"

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
