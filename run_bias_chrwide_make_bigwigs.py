import pandas as pd
import os

#data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_subsample_atac.csv",sep=",", names=["fold", "cell", "cell1", "model"])
#data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",", names=["fold", "cell", "model"])
data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",", names=["fold", "cell", "model"])
#print(data)
#data=pd.read_csv("k562_atac_uncorrected.csv",sep=",", names=["fold", "model"])
#data=pd.read_csv("k562_dnase_uncorrected.csv",sep=",", names=["fold", "model"])
print(data)

#output_dir=/mnt/lab_data2/anusri/print_analysis/make_bg_bigwig/$cellline/$dtype"_bias"
#bias_chrwide_make_bigwig.sh

for i,r in data.iterrows():

        if r["fold"] != "fold_0":
                continue

        if r["cell"] != "GM12878":
                continue

        chrombpnet_nb=r["model"]+"/chrombpnet_model/bias_model_scaled.h5"
        chrombpnet=r["model"]+"/chrombpnet_model/bias_model_scaled.h5"
        #chrombpnet=r["model"]
        #chrombpnet_nb=r["model"]
        #cellline=r["cell"]
        cellline="GM12878"
        outputf="/mnt/lab_data2/anusri/print_analysis/make_bg_bigwig/"+cellline+"/ATAC_bias_"+r["fold"].replace("_","")
        #"results/chrombpnet/auprc_curves/"+cellline+"/"+"/ATAC_uncorrected/"+cellline+"_"+r["fold"]
        #gpu="MIG-166d7783-762d-5f61-b31c-549eb4e0fba0"
        gpu="MIG-f80e9374-504a-571b-bac0-6fb00750db4c"
        #gpu="MIG-f80e9374-504a-571b-bac0-6fb00750db4c"
        #gpu="1"

        #print(chrombpnet_nb)
        if os.path.isfile(chrombpnet_nb):
                #ofile=outputf+"_wo_bias_predictions.h5"
                ofile=outputf+"_bias.bw"
                print(ofile)
                if not os.path.isfile(ofile):
                        command = "bash bias_chrwide_make_bigwig.sh "+chrombpnet_nb+" "+chrombpnet+" "+cellline+" "+r["fold"].replace("_","")+" ATAC "+gpu
                        print(command)
                        os.system(command)
                else:
                        print("Done !!!!!!!!!    "+ofile)

