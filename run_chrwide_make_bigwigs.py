import pandas as pd
import os

#data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_subsample_atac.csv",sep=",", names=["fold", "cell", "cell1", "model"])
#data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",", names=["fold", "cell", "model"])
#data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_atac.csv",sep=",", names=["fold", "cell", "model"])
#print(data)
data=pd.read_csv("k562_atac_uncorrected.csv",sep=",", names=["fold", "model"])
#data=pd.read_csv("k562_dnase_uncorrected.csv",sep=",", names=["fold", "model"])
print(data)



for i,r in data.iterrows():

        #if r["fold"] == "fold_0":
        #        continue

        #chrombpnet_nb=r["model"]+"/chrombpnet_model/chrombpnet_wo_bias.h5"
        #chrombpnet=r["model"]+"/chrombpnet_model/chrombpnet.h5"
        chrombpnet=r["model"]
        chrombpnet_nb=r["model"]
        #cellline=r["cell"]
        cellline="K562"
        outputf="results/chrombpnet/auprc_curves/"+cellline+"/"+"/ATAC_uncorrected/"+cellline+"_"+r["fold"]
        #gpu="MIG-166d7783-762d-5f61-b31c-549eb4e0fba0"
        #gpu="MIG-f80e9374-504a-571b-bac0-6fb00750db4c"
        gpu="MIG-f80e9374-504a-571b-bac0-6fb00750db4c"
        #gpu="3"

        #print(chrombpnet_nb)
        if os.path.isfile(chrombpnet_nb):
                #ofile=outputf+"_wo_bias_predictions.h5"
                ofile=outputf+"_wo_bias.bw"
                print(ofile)
                if not os.path.isfile(ofile):
                        command = "bash chr1wide_make_bigwigs.sh "+chrombpnet_nb+" "+chrombpnet+" "+cellline+" "+r["fold"].replace("_","")+" ATAC "+gpu
                        print(command)
                        os.system(command)
                else:
                        print("Done !!!!!!!!!    "+ofile)

