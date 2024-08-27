import pandas as pd
import os

#data = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_atac.csv",header=None)
data = pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",header=None)
cell_types=["K562"]
temp_r="/oak/stanford/groups/akundaje/projects/chromatin-atlas-2022/chrombpnet/folds/DNASE/"+cell_types[0]

NARROWPEAK_SCHEMA = ["chr", "start", "end", "1", "2", "3", "4", "5", "6", "summit"]
itype="profile"
gpu="0"
fold=4

for cell_type in cell_types:
        ndata = data[data[1]==cell_type].reset_index()
        for i,r in ndata.iterrows():
                print(i,r[2])

                #if i!=fold:
                #        continue         
                beds_path = os.path.join(r[2],"chrombpnet_model/interpret/full_"+cell_type+".interpreted_regions_"+itype+".bed")
                if os.path.exists(beds_path):
                        beds_path = os.path.join(r[2],"chrombpnet_model/interpret/full_"+cell_type+".interpreted_regions_"+itype+".bed")
                elif os.path.exists(os.path.join(r[2],"chrombpnet_model/interpret/merged."+cell_type+".interpreted_regions.bed")):
                        beds_path = os.path.join(r[2],"chrombpnet_model/interpret/merged."+cell_type+".interpreted_regions.bed")
                else:
                        beds_path = os.path.join(r[2],"interpret/merged."+cell_type+".interpreted_regions.bed")


                mdl_path = os.path.join(r[2],"chrombpnet_model/chrombpnet_wo_bias.h5")
                if not os.path.exists(mdl_path):
                        break

                mdl_path = os.path.join(r[2],"chrombpnet_model/")
                try:
                    if not os.path.isdir(os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/")):
                        os.mkdir(os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/"))
                    output_prefix=os.path.join(r[2],"chrombpnet_model/predictions_all_new_jan_2024/all_regions_preds")
                except:
                    if not os.path.isdir(os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/")):
                        os.mkdir(os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/"))
                    output_prefix=os.path.join(temp_r,"fold_"+str(i)+"_predictions_all_new_jan_2024/all_regions_preds")

                command = "bash script.sh "+ mdl_path+" "+beds_path+" "+output_prefix+" "+gpu
                print(command)
                os.system(command)



