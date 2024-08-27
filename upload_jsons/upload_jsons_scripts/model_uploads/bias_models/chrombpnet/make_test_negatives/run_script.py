import pandas as pd
import os

model_atac = pd.read_csv("../atac_bias_model_chrombpnet.csv",sep=",", header=None)


print(model_atac.head())

for i,r in model_atac.iterrows():

	print(r)
	if os.path.isfile(os.path.join(r[2], "train_test_regions_may_7_2024/nonpeaks.testset.bed.gz")):
		try:
			tdata = pd.read_csv(os.path.join(r[2], "train_test_regions_may_7_2024/nonpeaks.testset.bed.gz"))
			continue
		except:
			pass

	print(os.path.join(r[2], "train_test_regions_may_7_2024/nonpeaks.testset.bed.gz"))
	if r[0].split("_")[-1] == "0":
		command = "bash script.sh "+r[2]+" "+r[1]+" "+r[0]
	else:
		command = "bash script.sh "+r[2]+" "+r[1]+" "+r[0]+" "+"_"+str(r[0].split("_")[-1])

	print(command)
	os.system(command)	
