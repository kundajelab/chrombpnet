import pandas as pd
import os

data=pd.read_csv("logs/checkpoint/JAN_02_2023/model_dir_dnase.csv",sep=",", names=["fold", "cell", "model"])
#print(data)

for i,r in data.iterrows():
	if r["cell"] != "H1ESC":
		continue


	chrombpnet_nb=r["model"]+"/chrombpnet_model/chrombpnet_wo_bias.h5"
	cellline=r["cell"]
	model_dir=r["model"]
	gpu="5"

	#print(chrombpnet_nb)
	if os.path.isfile(chrombpnet_nb):
		ofile=model_dir+"/"+"background_interpret_new/"+cellline+".profile.bw"
		if not os.path.isfile(ofile):
			command = "bash make_bigwig_nonpeaks_all.sh "+chrombpnet_nb+" "+cellline+" "+model_dir+" "+gpu+" "+"DNASE_SE"
			print(command)
			os.system(command)
		else:
			print("Done !!!!!!!!!    "+ofile)
