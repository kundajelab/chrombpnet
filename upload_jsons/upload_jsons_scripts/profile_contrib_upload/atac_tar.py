import os
import json
import pandas as pd

#encids = ["K562", "HEPG2", "IMR90_new", "H1ESC_new", "GM12878_new"]
encids = ["IMR90", "H1ESC", "GM12878"]

encode_id = {"K562": "ENCSR868FGK",
"GM12878": "ENCSR637XSC",
"HEPG2": "ENCSR291GJU",
"IMR90": "ENCSR200OML",
"H1ESC": "GSE267154"} 
odir='atac/'

model_atac = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",", header=None)

def fetch_per_fold_profile(odir,model_path, encid, i, name):

		model_path_orig=model_path
		model_path="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/"+model_path.split("/")[-1]
		data_paths = []
		log_paths = []
		log_paths_opt = []
		
		odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/interpret_upload/fold_"+str(i)+"/"
		input_h5 = os.path.join(odir, name+"_profile_attribs_reformatted.h5")
		if os.path.isfile(input_h5):
			data_paths.append((input_h5, "seq_contrib.profile.fold_"+str(i)+".all_regions."+encid+".h5"))
		
		#model_path="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/ATAC_SE_04.27.2024//chrombpnet_model"
		
		# ATAC regions logs
		
		model_path=model_path+"/chrombpnet_model"
		input_log=model_path+"/interpret_dnase/full_"+name+".interpret.args.json"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.dnase_regions.fold_"+str(i)+"."+encid+".args.json"))
		else:
		
			print(input_log)
		input_log=model_path+"/interpret_dnase/full_"+name+".interpet.log"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.dnase_regions.fold_"+str(i)+"."+encid+".log"))
		else:
			print(input_log)

		input_log=model_path+"/interpret_dnase/ATAC_peaks_full.profile.interpret.log.e"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.dnase_regions.fold_"+str(i)+"."+encid+".stderr.txt"))
		else:
			print(input_log)

		input_log=model_path+"/interpret_dnase/ATAC_peaks_full.profile.interpret.log.o"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.dnase_regions.fold_"+str(i)+"."+encid+".stdout.txt"))
		else:
			print(input_log)

		# atac regions logs
		
		input_log=model_path_orig+"/interpret/merged."+name+".interpret.args.json"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.atac_regions.fold_"+str(i)+"."+encid+".args.json"))
		else:
			print(input_log)

		input_log=model_path_orig+"/interpret/merged."+name+".interpet.log"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.atac_regions.fold_"+str(i)+"."+encid+".log"))
		else:
			print(input_log)
			
		# atac regions logs


		input_log=model_path+"/interpret/full_"+name+".interpret.args.json"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.atac_regions.fold_"+str(i)+"."+encid+".args.json"))
		else:
			print(input_log)

		input_log=model_path+"/interpret/full_"+name+".interpet.log"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.atac_regions.fold_"+str(i)+"."+encid+".log"))
		else:
			print(input_log)
			
		input_log=model_path+"/interpret/full.profile.interpret.log.e"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.atac_regions.fold_"+str(i)+"."+encid+".stderr.txt"))
		else:
			print(input_log)

		input_log=model_path+"/interpret/full.profile.interpret.log.o"
		if os.path.isfile(input_log):
			log_paths.append((input_log, "logs.seq_contrib.profile.atac_regions.fold_"+str(i)+"."+encid+".stdout.txt"))
		else:
			print(input_log)
			
                          
		return data_paths, log_paths, log_paths_opt
        
def fetch_profile_tar(encid, args_json, model_paths, name):
		success = False
		args_json["profile sequence contribution scores tar"] = {}
		readme_file = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/dummy/chrombpnet_test/READMEs/profile.deepshap.README"
		assert(os.path.isfile(readme_file))
		args_json["profile sequence contribution scores tar"]["file.paths"] = [(readme_file, "README.md")]
		args_json["profile sequence contribution scores tar"]["logs.seq_contrib.profile."+encid] = {"file.paths": []}
		
		## full h5 path
		
		odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/interpret_upload/average_preds/"
		
		input_h5 = os.path.join(odir, name+"_profile_attribs_reformatted.h5")
		if os.path.isfile(input_h5):
				args_json["profile sequence contribution scores tar"]["file.paths"].append((input_h5,"seq_contrib.profile.fold_mean.all_regions."+encid+".h5"))               
		else:
				success = False
				return success, args_json
		
		## modisoc h5 path
		
		modisco_input = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/merge_folds_new_may_05_24/in_peaks.profile_scores_new_compressed.h5"
		if os.path.isfile(modisco_input):
				args_json["profile sequence contribution scores tar"]["file.paths"].append((modisco_input,"seq_contrib.profile.fold_mean.modisco_input."+encid+".h5"))               
		else:
				success = False
				return success, args_json
		
		# log files 
		
		
		input_file=model_paths[0]+"/chrombpnet_model/interpret_all/full_"+name+".interpreted_regions_profile.bed"
		newf="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/interpret_upload/average_preds/per_folds.inputs.bed.gz"
		input_bed = pd.read_csv(input_file, compression='gzip', sep='\t', header=None) 
		if os.path.isfile(input_file):
			if not os.path.isfile(newf):
				input_bed.to_csv(newf, sep='\t', header=False, index=False, compression='gzip')
			args_json["profile sequence contribution scores tar"]["logs.seq_contrib.profile."+encid]["file.paths"].append((newf,"logs.seq_contrib.profile.input_regions.per_fold."+encid+".bed.gz"))              
            
			newf="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/interpret_upload/average_preds/mean_folds.inputs.bed.gz"
			if not os.path.isfile(newf):
				input_bed = input_bed[~(input_bed[0]=="chrM")]
				input_bed.to_csv(newf, sep='\t', header=False, index=False, compression='gzip')
			args_json["profile sequence contribution scores tar"]["logs.seq_contrib.profile."+encid]["file.paths"].append((newf,"logs.seq_contrib.profile.input_regions.fold_mean."+encid+".bed.gz"))   
		
		input_file="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/merge_folds_new_may_05_24/in_peaks.profile.interpreted_regions.bed"
		newf="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/interpret_upload/average_preds/modisco.inputs.bed.gz"
		input_bed = pd.read_csv(input_file, sep='\t', header=None) 
		if os.path.isfile(input_file):
			if not os.path.isfile(newf):
				input_bed.to_csv(newf, sep='\t', header=False, index=False, compression='gzip')
			args_json["profile sequence contribution scores tar"]["logs.seq_contrib.profile."+encid]["file.paths"].append((newf,"logs.seq_contrib.profile.input_regions.modisco."+encid+".bed.gz"))              
		
		odir="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+name+"/interpret_upload/average_preds/"
		
		input_log = os.path.join(odir, "reformat.log.e")
		if os.path.isfile(input_log):
				args_json["profile sequence contribution scores tar"]["logs.seq_contrib.profile."+encid]["file.paths"].append((input_log, "logs.seq_contrib.profile.fold_mean.reformat."+encid+".stderr.txt"))
		
		input_log = os.path.join(odir, "reformat.log.o")
		if os.path.isfile(input_log):
				args_json["profile sequence contribution scores tar"]["logs.seq_contrib.profile."+encid]["file.paths"].append((input_log, "logs.seq_contrib.profile.fold_mean.reformat."+encid+".stdout.txt"))
			   
		assert(len(args_json["profile sequence contribution scores tar"]["logs.seq_contrib.profile."+encid]["file.paths"])==5) 
						
		for i in range(5):
				data_paths, log_paths, log_paths_opt = fetch_per_fold_profile(odir,model_paths[i], encid, i, name)
		
				if data_paths is None:
						success = False
						return success, args_json
						
				args_json["profile sequence contribution scores tar"]["fold_"+str(i)] = {}
				args_json["profile sequence contribution scores tar"]["fold_"+str(i)]["file.paths"] = data_paths
				args_json["profile sequence contribution scores tar"]["fold_"+str(i)]["logs.seq_contrib.profile.fold_"+str(i)+"."+encid] = {"file.paths": log_paths+log_paths_opt}
				assert(len(data_paths) == 1)
				print(len(log_paths))
				assert(len(log_paths) >= 4)
								
		success=True
		return success, args_json

for encid in encids:
		print(encid)
		
		
		ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores_new_compressed.stats"
		if os.path.isfile(ofile):
				profile_bw = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores_new_compressed.bw"
		else:
				profile_bw = None
				print(ofile)
				
		ofile = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores_new_compressed.stats"
		if os.path.isfile(ofile):
				profile_bw = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/ATAC/"+encid+"/interpret_upload/average_preds/"+encid+"_folds_merged.profile_scores_new_compressed.bw"
		else:   
				profile_bw = None               
				print(ofile)
				continue
		
		assert(os.path.isfile(profile_bw)==True) 
		assert(os.path.isfile(profile_bw)==True)        
		
		model_paths = model_atac[model_atac[1]==encid.replace("_new","")][2].values
		print(model_paths)
		args_json = {}
		args_json["experiment"] = encode_id[encid]
		
		
		success, args_json = fetch_profile_tar(encode_id[encid], args_json, model_paths, encid)
		if not success:
				print("ERR profile tar")
				continue                

		if not os.path.isfile(odir+encode_id[encid]+".json"):
				f = open(odir+encode_id[encid]+".json", "w")
				json.dump(args_json, f, indent=4)
				f.close()
