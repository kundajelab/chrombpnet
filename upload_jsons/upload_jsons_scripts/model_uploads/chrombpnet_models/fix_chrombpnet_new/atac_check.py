import os
import json
import pandas as pd

odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/ATAC/"
bw_odir = "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/full_deepshaps/bigwigs/ATAC/"
#output_dir =  "/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022-uploads/jsons/ATAC/stage1/jul_26_2023/"


models_path = ["chrombpnet_model_feb15", "chrombpnet_model_feb15_fold_1", "chrombpnet_model_feb15_fold_2", "chrombpnet_model_feb15_fold_3", "chrombpnet_model_feb15_fold_4"]
encids = open("../data/atac_passed.txt").readlines()
encids = [line.strip() for line in encids]

list_of_models = []

ignore_list = ["ENCSR637XSC", "ENCSR291GJU", "ENCSR200OML", "ENCSR868FGK"]

for encid in encids:

	if encid in ignore_list:
		continue
			
	for i in range(5):
		model_dir=models_path[i]
		checks_file = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/check_passed.txt")
		if os.path.isfile(checks_file):
				continue
		
		cm_model = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet_new.tar")
		cm_model_new = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet_new/")
		cm_model_new_v2 = os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet/")
		#print(os.path.isfile(os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar")))
		#print(os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar"))
		if not os.path.isfile(os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/chrombpnet.tar")):
			if os.path.isfile(cm_model):

				#command = "cp -r "+cm_model_new+" "+cm_model_new_v2
				#print(command)
				#os.system(command)
				list_of_models.append([encid, model_dir])
				#output_dir=os.path.join(odir, encid + "/" + model_dir + "/new_chrombpnet_model/")
				#file_path="chrombpnet"
				#command = "cd "+output_dir+" && tar -cf "+file_path+".tar "+file_path+"/"
				#print(command)
				#os.system(command)
			
			else:

				print("no tar file! "+encid, model_dir)	
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
df_models.to_csv("new_dir_copied_atac.txt", sep=",", header=False, index=False)
