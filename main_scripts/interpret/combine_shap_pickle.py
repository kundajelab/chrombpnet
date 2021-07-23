import argparse
import os
import pickle as pkl

parser=argparse.ArgumentParser(description="view model arch")
parser.add_argument("--source")
parser.add_argument("--target")
parser.add_argument("--type",type=str)
args=parser.parse_args()


if args.type=="20k":
        f1 = os.path.join(args.source,"scores.xaa.fold0.deepSHAP")
        f2 = os.path.join(args.source,"scores.xab.fold0.deepSHAP")
        f3 = os.path.join(args.source,"scores.xac.fold0.deepSHAP")
        fileso = [f1, f2, f3]
else:
        f1 = os.path.join(args.source,"scores.xaa.fold0.deepSHAP")
        f2 = os.path.join(args.source,"scores.xab.fold0.deepSHAP")
        f3 = os.path.join(args.source,"scores.xac.fold0.deepSHAP")
        f4 = os.path.join(args.source,"scores.xad.fold0.deepSHAP")
        f5 = os.path.join(args.source,"scores.xae.fold0.deepSHAP")
        f6 = os.path.join(args.source,"scores.xaf.fold0.deepSHAP")
        f7 = os.path.join(args.source,"scores.xag.fold0.deepSHAP")
        f8 = os.path.join(args.source,"scores.xah.fold0.deepSHAP")
        f9 = os.path.join(args.source,"scores.xai.fold0.deepSHAP")
        f10 = os.path.join(args.source,"scores.xaj.fold0.deepSHAP")
        fileso = [f1, f2, f3, f4, f5, f6, f7, f8, f9, f10]

data1 = pkl.load(open(f1, "rb"))
data3 = {}
data3["seq"] = {}
data3["count_shap"] = {}
data3["profile_shap"] = {}

for fin in fileso:
	data2 = pkl.load(open(fin, "rb"))
	for key in data2["seq"]:
		if key not in data3["seq"]:
			data3["seq"][key] = data2["seq"][key]
			data3["count_shap"][key] = data2["count_shap"][key]
			data3["profile_shap"][key] = data2["profile_shap"][key]
		else:
			print(key)

print(len(data3["seq"].keys()))

if args.type == "20k":
	f3 = os.path.join(args.target,"gm12878.20K.dnase.fold0.deepSHAP")
else:
	f3 = os.path.join(args.target,"gm12878.all.dnase.fold0.deepSHAP")

of = open(f3,"wb")
pkl.dump(data3,of)
