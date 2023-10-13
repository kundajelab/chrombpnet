import os
import json
import pandas as pd

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/DNASE/"
bw_odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/full_deepshaps/bigwigs/DNASE/"
#output_dir =  "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/DNASE/stage1/jul_26_2023/"
output_dir="dnase_production_uploads/"

tissue_encids = open("../data/tissue_passed.txt").readlines()
tissue_encids = [line.strip() for line in tissue_encids]

primary_encids = open("../data/primary_passed.txt").readlines()
primary_encids = [line.strip() for line in primary_encids]

celline_encids = open("../data/cellline_passed.txt").readlines()
celline_encids = [line.strip() for line in celline_encids]

invitro_encids = open("../data/invitro_passed.txt").readlines()
invitro_encids = [line.strip() for line in invitro_encids]

primary_models_path = ["chrombppnet_model_encsr283tme_bias", "chrombppnet_model_encsr283tme_bias_fold_1", "chrombppnet_model_encsr283tme_bias_fold_2", "chrombppnet_model_encsr283tme_bias_fold_3", "chrombppnet_model_encsr283tme_bias_fold_4"]
celline_models_path = ["chrombpnet_model_feb15_fold_0", "chrombpnet_model_feb15_fold_1", "chrombpnet_model_feb15_fold_2", "chrombpnet_model_feb15_fold_3", "chrombpnet_model_feb15_fold_4"]
tissue_models_path = ["chrombpnet_model_encsr880cub_bias","chrombppnet_model_encsr880cub_bias_fold_1","chrombppnet_model_encsr880cub_bias_fold_2","chrombppnet_model_encsr880cub_bias_fold_3","chrombppnet_model_encsr880cub_bias_fold_4"]
invitro_models_path = ["chrombpnet_model_encsr146kfx_bias", "chrombpnet_model_encsr146kfx_bias_fold_1", "chrombpnet_model_encsr146kfx_bias_fold_2", "chrombpnet_model_encsr146kfx_bias_fold_3", "chrombpnet_model_encsr146kfx_bias_fold_4"]

encids = tissue_encids + primary_encids + celline_encids + invitro_encids

list_of_models = []

for encid in encids:

	ignore_list = ["ENCSR000EMT", "ENCSR149XIL", "ENCSR477RTP", "ENCSR000EOT", "ENCSR000EMU"]

	if encid  in ignore_list:
		continue
	
	if encid in primary_encids:
		models_path = primary_models_path
		bias_encid="ENCSR283TME"
		#print("primary")
	elif encid in tissue_encids:
		models_path = tissue_models_path
		bias_encid="ENCSR880CUB"
		#print("tissue")
	elif encid in invitro_encids:
		models_path = invitro_models_path
		bias_encid="ENCSR146KFX"
		#print("invitro")
	elif encid in celline_encids:
		models_path = celline_models_path
		bias_encid="ENCSR149XIL"
		#print("celline")
	else:
		print(encid)
		print("type not found")
		continue
			
	for i in range(5):
		model_dir=models_path[i]
		checks_file = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/check_passed.txt")
		if os.path.isfile(checks_file):
				continue
		
		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet_new.tar")
		cm_model_new = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet_new/")
		cm_model_new_v2 = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet/")
		print(os.path.isfile(os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar")))
		print(os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar"))
		if not os.path.isfile(os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar")):
			if os.path.isfile(cm_model):

				command = "cp -r "+cm_model_new+" "+cm_model_new_v2
				print(command)
				os.system(command)
				list_of_models.append([encid, model_dir])
				output_dir=os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/")
				file_path="chrombpnet"
				command = "cd "+output_dir+" && tar -cf "+file_path+".tar "+file_path+"/"
				print(command)
				os.system(command)
			
			else:

				print("no tar file! "+encid+model_dir)	
		else:

			#command = "rm -r "+os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet/*")
			#os.system(command)
			#print(command)

			#command = "rm "+os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar")
			#os.system(command)
			#print(command)

			if os.path.isfile(cm_model):
				list_of_models.append([encid, model_dir])
			
df_models = pd.DataFrame(list_of_models)
df_models.to_csv("new_dir_copied.txt", sep=",", header=False, index=False)
