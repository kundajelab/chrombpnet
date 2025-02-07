import os
import json
import pandas as pd

#names = ["IMR90_new", "H1ESC_new", "GM12878_new"]
names = ["K562", "HEPG2"]


#model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/v1/model_dir_dnase_v2_interpret.csv",sep=",", header=None)
model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/v1/model_dir_dnase_v2.csv",sep=",", header=None)

encode_id = {"HEPG2": "ENCSR149XIL",
        "K562": "ENCSR000EOT",
        "IMR90_new": "ENCSR477RTP",
        "GM12878_new": "ENCSR000EMT",
        "H1ESC_new": "ENCSR000EMU"}


encode_id_dnase = {
"GM12878_new": "ENCSR000EMT",
"IMR90_new": "ENCSR477RTP",
"H1ESC_new": "ENCSR000EMU"}  

outdir='dnase_tar/'

def fetch_per_fold_preds(odir,model_path, encid, i, name):

	data_paths = []
	log_paths = []
	log_paths_opt = []
	
	odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/preds_upload/fold_"+str(i)+"/"
	input_h5 = os.path.join(odir, name+"_wo_bias_all_with_ccre_predictions.h5")
	if os.path.isfile(input_h5):
		data_paths.append((input_h5, "pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".h5"))

	input_log = os.path.join(odir, "pred.counts.log.e")
	#print(input_log)
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stderr.txt"))

	input_log = os.path.join(odir, "pred.counts.log.o")
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stdout.txt"))
		
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/"+model_path.split("/")[-2]+"/chrombpnet_model/preds_atac/pred.counts.log.o"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stdout_v1.txt"))

	#print(input_log)
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/"+model_path.split("/")[-2]+"/chrombpnet_model/preds_atac/pred.counts.log.e"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stderr_v1.txt"))
		
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/"+model_path.split("/")[-2]+"/chrombpnet_model/preds_dnase/pred.counts.log.o"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stdout_v2.txt"))
		
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/"+model_path.split("/")[-2]+"/chrombpnet_model/preds_dnase/pred.counts.log.e"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stderr_v2.txt"))
	
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/"+model_path.split("/")[-2]+"/chrombpnet_model/preds_ccre/pred.counts.log.o"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stdout_v3.txt"))
		
	input_log="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/"+model_path.split("/")[-2]+"/chrombpnet_model/preds_ccre/pred.counts.log.e"
	if os.path.isfile(input_log):
		log_paths.append((input_log, "logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid+".stderr_v3.txt"))
		
						
	return data_paths, log_paths, log_paths_opt

def fetch_pred_tar(encid, args_json, model_paths, name):
	success = False
	args_json["bias-corrected predicted signal profile tar"] = {}
	readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/dummy/chrombpnet_test/READMEs/bc.predicted.README"
	assert(os.path.isfile(readme_file))
	args_json["bias-corrected predicted signal profile tar"]["file.paths"] = [(readme_file, "README.md")]
	args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid] = {"file.paths": []}

	odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/preds_upload/average_preds_with_ccre_vf/"

	input_h5 = os.path.join(odir, name+".mean_preds_wo_bias_predictions.h5")
	if os.path.isfile(input_h5):
		args_json["bias-corrected predicted signal profile tar"]["file.paths"].append((input_h5,"pred.chrombpnet_nobias.fold_mean."+encid+".h5"))		
	else:
		success = False
		return success, args_json

	if name in ["IMR90_new", "GM12878_new", "H1ESC_new"]:
		bed1 = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"+encode_id_dnase[name]+"/preprocessing/downloads/peaks.bed.gz", sep='\t', header=None)
	else:	
		bed1 = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/peaks_no_blacklist.bed", sep='\t', header=None)

	
	bed2 = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name.replace("_new","")+"/peaks_no_blacklist.bed", sep='\t', header=None)

	bed3 = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/ccres/ccres_fill.bed.gz", sep='\t', header=None)

	print(bed1.shape)
	print(bed2.shape)
	bedf = pd.concat([bed2, bed1, bed3])
	print(bedf.shape)
	
# 	input_bed = os.path.join(odir, "input.regions.bed.gz")
# 	if os.path.isfile(input_bed):
# 		args_json["bias-corrected predicted signal profile tar"]["file.paths"].append((input_bed,"input_regions.pred.chrombpnet_nobias."+encid+".bed.gz"))		
# 	else:
# 		bedf.to_csv(input_bed, sep='\t', header=False, index=False, compression='gzip')

	input_bed = os.path.join(odir, "input.regions.bed.gz")
	if os.path.isfile(input_bed):
		bedf.to_csv(input_bed, sep='\t', header=False, index=False, compression='gzip')
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_bed,"logs.pred.chrombpnet_nobias.input_regions.per_fold."+encid+".bed.gz"))		
	else:
		bedf.to_csv(input_bed, sep='\t', header=False, index=False, compression='gzip')
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_bed,"logs.pred.chrombpnet_nobias.input_regions.per_fold."+encid+".bed.gz"))		

	bedn = pd.read_csv("/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/preds_upload/fold_0/"+name+"_wo_bias_all_with_ccre_regions.bed", sep='\t', header=None)

	input_bed = os.path.join(odir, "filtered.regions.bed.gz")
	if os.path.isfile(input_bed):
		bedn.to_csv(input_bed, sep='\t', header=False, index=False, compression='gzip')

		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_bed,"logs.pred.chrombpnet_nobias.input_regions.fold_mean."+encid+".bed.gz"))		
	else:
		bedn.to_csv(input_bed, sep='\t', header=False, index=False, compression='gzip')
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_bed,"logs.pred.chrombpnet_nobias.input_regions.fold_mean."+encid+".bed.gz"))		


	input_log = os.path.join(odir, "merge.preds.log.e")
	if os.path.isfile(input_log):
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_log, "logs.pred.chrombpnet_nobias.fold_mean."+encid+".stderr.txt"))

	input_log = os.path.join(odir, "merge.preds.log.o")
	if os.path.isfile(input_log):
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_log, "logs.pred.chrombpnet_nobias.fold_mean."+encid+".stdout.txt"))
		
	input_log = os.path.join(odir, name+".mean_preds_wo_bias.stat")
	if os.path.isfile(input_log):
		args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"].append((input_log, "logs.pred.chrombpnet_nobias.fold_mean."+encid+".stat"))

	assert(len(args_json["bias-corrected predicted signal profile tar"]["logs.pred.chrombpnet_nobias.fold_mean."+encid]["file.paths"])==5)
		
		
	for i in range(5):
		print(i)
		data_paths, log_paths, log_paths_opt = fetch_per_fold_preds(odir,model_paths[i], encid, i, name)

		if data_paths is None:
			success = False
			return success, args_json
			
		args_json["bias-corrected predicted signal profile tar"]["fold_"+str(i)] = {}
		args_json["bias-corrected predicted signal profile tar"]["fold_"+str(i)]["file.paths"] = data_paths
		args_json["bias-corrected predicted signal profile tar"]["fold_"+str(i)]["logs.pred.chrombpnet_nobias.fold_"+str(i)+"."+encid] = {"file.paths": log_paths+log_paths_opt}
		assert(len(data_paths) == 1)
		print(len(log_paths))
		assert(len(log_paths) == 8)
		#assert(len(log_paths) >= 2)
				
	success=True
	return success, args_json


for name in names:
	ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+name+"/preds_upload/average_preds_with_ccre_vf/"+name+".mean_preds_wo_bias.stat"
	if os.path.isfile(ofile):
		args_json = {}
		encid=encode_id[name]
		args_json["experiment"] = encid
		model_paths = model_atac[model_atac[1]==name.replace("_new","")][2].values
		print(model_paths)
		success, args_json = fetch_pred_tar(encid, args_json, model_paths, name)
		if not success:
			print("ERR preds tar")
			continue		
			
		if not os.path.isfile(outdir+encid+".json"):
			f = open(outdir+encode_id[name]+".json", "w")
			json.dump(args_json, f, indent=4)
			f.close()

