import argparse
import json
import numpy as np
from scipy import stats

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--files", type=str)
args = parser.parse_args()


files1 = args.files.strip().split(",")
print(files1)
pr=[]
sr=[]
mse=[]
jsd=[]

for file in files1:
	data = json.load(open(file+"/chrombpnet_model/chrombpnet_metrics.json"))
	pr.append(data["counts_metrics"]["peaks"]["pearsonr"])
	sr.append(data["counts_metrics"]["peaks"]["spearmanr"])
	mse.append(data["counts_metrics"]["peaks"]["mse"])
	jsd.append(data["profile_metrics"]["peaks"]["median_norm_jsd"])

pr8=[]
sr8=[]
mse8=[]
jsd8=[]


files2 = ["HEPG2_05.27.2022_1234_8_2114_0_only_pos","HEPG2_05.25.2022_1234_8_2114_1_only_pos", "HEPG2_05.25.2022_1234_8_2114_2_only_pos", "HEPG2_05.27.2022_1234_8_2114_3_only_pos", "HEPG2_05.27.2022_1234_8_2114_4_only_pos"]
for file in files2:
	data = json.load(open(file+"/chrombpnet_model/chrombpnet_metrics.json"))
	pr8.append(data["counts_metrics"]["peaks"]["pearsonr"])
	sr8.append(data["counts_metrics"]["peaks"]["spearmanr"])
	mse8.append(data["counts_metrics"]["peaks"]["mse"])
	jsd8.append(data["profile_metrics"]["peaks"]["median_norm_jsd"])



def get_result(st):
	a = ""
	a+=str(np.round(st[0],4))+", "
	a+=str(np.round(st[1],4))
	return a

print("Wilcoxon")
print("pearsonr",get_result(stats.wilcoxon(pr,pr8)))
print("spearmanr", get_result(stats.wilcoxon(sr,sr8)))
print("mse",get_result(stats.wilcoxon(mse,mse8)))
print("jsd",get_result(stats.wilcoxon(jsd,jsd8)))
print(jsd)
print(jsd8)
print(pr)
print(pr8)
print(mse)
print(mse8)



print("paired t-test")
print("pearsonr",get_result(stats.ttest_rel(pr,pr8)))
print("spearmanr",get_result(stats.ttest_rel(sr,sr8)))
print("mse",get_result(stats.ttest_rel(mse,mse8)))
print("jsd",get_result(stats.ttest_rel(jsd,jsd8)))


