import os
import json
import pandas as pd

names = ["K562", "HEPG2", "IMR90", "H1ESC", "GM12878"]

encode_id = {"K562": "ENCSR868FGK",
"GM12878": "ENCSR637XSC",
"HEPG2": "ENCSR291GJU",
"IMR90": "ENCSR200OML",
"H1ESC": "GSE267154"}  

encode_id_dnase = {
"GM12878": "ENCSR000EMT",
"IMR90": "ENCSR477RTP",
"H1ESC": "ENCSR000EMU"}  

outdir='atac_tar/'

def fetch_per_fold_preds(odir,model_path, encid, i, name):

	data_paths = []
	log_paths = []
	log_paths_opt = []
	
	odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/preds_upload/fold_"+str(i)+"/"
	input_h5 = os.path.join(odir, name+"_w_bias_all_predictions.h5")
	data_paths.append((input_h5, "pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".h5"))

	input_log = os.path.join(odir, "pred.counts.log.e")
	#print(input_log)
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stderr.txt"))

	input_log = os.path.join(odir, "pred.counts.log.o")
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stdout.txt"))
		
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/"+model_path.split("/")[-1]+"/chrombpnet_model/preds_atac/pred.counts.log.o"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stdout_v1.txt"))

	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/"+model_path.split("/")[-1]+"/chrombpnet_model/preds_atac/pred.counts.log.e"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stderr_v1.txt"))
		
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/"+model_path.split("/")[-1]+"/chrombpnet_model/preds_dnase/pred.counts.log.o"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stdout_v2.txt"))
		
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/"+model_path.split("/")[-1]+"/chrombpnet_model/preds_dnase/pred.counts.log.e"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stderr_v2.txt"))
					
	return data_paths, log_paths, log_paths_opt

def fetch_pred_tar(encid, args_json, model_paths, name):
	success = False
	args_json["bias-corrected predicted signal profile tar"] = {}
	readme_file = "READMEs/bc.predicted.README"
	assert(os.path.isfile(readme_file))
	args_json["bias-corrected predicted signal profile tar"]["file.paths"] = [(readme_file, "README.md")]
	args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid] = {"file.paths": []}

	odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/preds_upload/average_preds/"

	input_h5 = os.path.join(odir, name+".mean_preds_w_bias_predictions.h5")
	if os.path.isfile(input_h5):
		args_json["bias-corrected predicted signal profile tar"]["file.paths"].append((input_h5,"pred.chrombpnet_nobias.fold_mean."+encid+".h5"))		
	else:
		success = False
		return success, args_json

	if name in ["IMR90", "GM12878", "H1ESC"]:
		bed1 = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"+encode_id_dnase[name]+"/preprocessing/downloads/peaks.bed.gz", sep='\t', header=None)
	else:	
		bed1 = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/peaks_no_blacklist.bed", sep='\t', header=None)

	
	bed2 = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/peaks_no_blacklist.bed", sep='\t', header=None)

	print(bed1.shape)
	print(bed2.shape)
	bedf = pd.concat([bed1, bed2])
	print(bedf.shape)
	
	input_bed = os.path.join(odir, "input.regions.bed.gz")
	if os.path.isfile(input_bed):
		args_json["bias-corrected predicted signal profile tar"]["file.paths"].append((input_bed,"input_regions.pred.chrombpnet_nobias."+encid+".bed.gz"))		
	else:
		bedf.to_csv(input_bed, sep='\t', header=False, index=False, compression='gzip')


	input_log = os.path.join(odir, "merge.preds.log.e")
	if os.path.isfile(input_log):
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_log, "logs.pred.chrombpnet_nobias.fold_mean."+encid+".stderr.txt"))

	input_log = os.path.join(odir, "merge.preds.log.o")
	if os.path.isfile(input_log):
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_log, "logs.pred.chrombpnet_nobias.fold_mean."+encid+".stdout.txt"))
		
			
	for i in range(5):
		data_paths, log_paths, log_paths_opt = fetch_per_fold_preds(odir,model_paths[i], encid, i, name)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["bias-corrected predicted signal profile tar"]["fold_"+str(i)] = {}
		args_json["bias-corrected predicted signal profile tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["bias-corrected predicted signal profile tar"]["fold_"+str(i)]["logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid] = {"file.paths": log_paths+log_paths_opt}
		assert(len(data_paths) == 1)
		print(len(log_paths))
		assert(len(log_paths) == 6)
				
	success=True
	return success, args_json

model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",", header=None)

for name in names:
	ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/preds_upload/average_preds/"+name+".mean_preds_w_bias.stat"
	if os.path.isfile(ofile):
		args_json = {}
		encid=encode_id[name]
		args_json["experiment"] = encid
		model_paths = model_atac[model_atac[1]==name][2].values
		print(model_paths)
		success, args_json = fetch_pred_tar(encid, args_json, model_paths, name)
		if not success:
			print("ERR preds tar")
			continue		
			
		if not os.path.isfile(outdir+encid+"_w_bias.json"):
			f = open(outdir+encode_id[name]+".json", "w")
			json.dump(args_json, f, indent=4)
			f.close()

