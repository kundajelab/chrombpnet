import pandas as pd
import os

data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_subsample_atac.csv",sep=",", names=["fold", "cell", "cell1", "model"])
#data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",", names=["fold", "cell", "model"])
#print(data)


for i,r in data.iterrows():
        #if r["cell"] != "GM12878":
        #        continue
        #print(r)
        if r["cell"] == "572M":
                continue

        if r["fold"] != "fold_0":
                continue

        chrombpnet_nb=r["model"]+"/chrombpnet_model/chrombpnet_wo_bias.h5"
        chrombpnet=r["model"]+"/chrombpnet_model/chrombpnet.h5"
        cellline=r["cell"]
        outputf="results/chrombpnet/auprc_curves/"+cellline+"/"+cellline
        gpu="MIG-166d7783-762d-5f61-b31c-549eb4e0fba0"

        #print(chrombpnet_nb)
        if os.path.isfile(chrombpnet_nb):
                #ofile=outputf+"_wo_bias_predictions.h5"
                ofile=outputf+"_wo_bias.bw"
                print(ofile)
                if not os.path.isfile(ofile):
                        command = "bash chr1wide_make_bigwigs.sh "+chrombpnet_nb+" "+chrombpnet+" "+cellline+" "+gpu
                        print(command)
                        os.system(command)
                else:
                        print("Done !!!!!!!!!    "+ofile)

