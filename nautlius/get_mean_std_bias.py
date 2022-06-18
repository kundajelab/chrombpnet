import argparse
import json
import numpy as np

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
	data = json.load(open(file+"/bias_model/bias_metrics.json"))
	pr.append(data["counts_metrics"]["nonpeaks"]["pearsonr"])
	sr.append(data["counts_metrics"]["nonpeaks"]["spearmanr"])
	mse.append(data["counts_metrics"]["nonpeaks"]["mse"])
	jsd.append(data["profile_metrics"]["nonpeaks"]["median_norm_jsd"])

print("pearsonrr")
print(np.round(np.mean(pr),2))
print(np.round(np.std(pr),2))

print("spearmanr")
print(np.round(np.mean(sr),2))
print(np.round(np.std(sr),2))


print("mse")
print(np.round(np.mean(mse),2))
print(np.round(np.std(mse),2))

print("jsd")
print(np.round(np.mean(jsd),2))
print(np.round(np.std(jsd),2))

pr=[]
sr=[]
mse=[]
jsd=[]



for file in files1:
	data = json.load(open(file+"/chrombpnet_model/bias_metrics.json"))
	pr.append(data["counts_metrics"]["peaks"]["pearsonr"])
	sr.append(data["counts_metrics"]["peaks"]["spearmanr"])
	mse.append(data["counts_metrics"]["peaks"]["mse"])
	jsd.append(data["profile_metrics"]["peaks"]["median_norm_jsd"])


print(pr)
print("pearsonrr")
print(np.round(np.mean(pr),2))
print(np.round(np.std(pr),2))

print("spearmanr")
print(np.round(np.mean(sr),2))
print(np.round(np.std(sr),2))


print("mse")
print(np.round(np.mean(mse),2))
print(np.round(np.std(mse),2))

print("jsd")
print(np.round(np.mean(jsd),2))
print(np.round(np.std(jsd),2))



