import json
import argparse
from scipy import stats 
import numpy as np

data = {"GM12878": {"invivo":"GM12878_03.01.2022_bias_128_4_1234_0.4_fold_0" , "transfer": "GM12878_06.08.2022_1234_8_2114_0_hepg2_transfer_bias"},
	"K562": {"invivo":"K562_02.17.2022_bias_128_4_1234_0.5_fold_0" , "transfer":"K562_06.08.2022_1234_8_2114_0_hepg2_transfer_bias" },
	"IMR90": {"invivo":"IMR90_04.09.2022_bias_128_4_1234_0.4_fold_0" , "transfer":"IMR90_06.08.2022_1234_8_2114_0_hepg2_transfer_bias" },
	"H1ESC":{"invivo":"H1ESC_05.09.2022_bias_128_4_1234_0.8_fold_0" , "transfer":"H1ESC_06.08.2022_1234_8_2114_0_hepg2_transfer_bias" }}

parser = argparse.ArgumentParser()
parser.add_argument("--metric",type=str)
parser.add_argument("--peaks_nonpeaks",type=str)
args=parser.parse_args()

invivo = []
transfer = []

for cellline in data:
	ipath = cellline+"/"+data[cellline]["invivo"]+"/chrombpnet_model/chrombpnet_wo_bias_metrics.json"
	tpath = cellline+"/"+data[cellline]["transfer"]+"/chrombpnet_model/chrombpnet_wo_bias_metrics.json"
	idata = json.load(open(ipath))
	tdata = json.load(open(tpath))

	print(cellline)
	invivo.append(idata["counts_metrics"][args.peaks_nonpeaks][args.metric])
	transfer.append(tdata["counts_metrics"][args.peaks_nonpeaks][args.metric])

invivo1 = [np.round(val,2) for val in invivo]
transfer1 = [np.round(val,2) for val in transfer]


print(invivo1)
print(transfer1)

print(stats.ttest_rel(invivo,transfer))


invivo = []
transfer = []

for cellline in data:
	ipath = cellline+"/"+data[cellline]["invivo"]+"/chrombpnet_model/chrombpnet_wo_bias_metrics.json"
	tpath = cellline+"/"+data[cellline]["transfer"]+"/chrombpnet_model/chrombpnet_wo_bias_metrics.json"
	idata = json.load(open(ipath))
	tdata = json.load(open(tpath))

	print(cellline)
	invivo.append(idata["profile_metrics"][args.peaks_nonpeaks]["median_jsd"])
	transfer.append(tdata["profile_metrics"][args.peaks_nonpeaks]["median_jsd"])

invivo1 = [np.round(val,2) for val in invivo]
transfer1 = [np.round(val,2) for val in transfer]



print(invivo1)
print(transfer1)

print(stats.ttest_rel(invivo,transfer))


invivo = []
transfer = []

for cellline in data:
	ipath = cellline+"/"+data[cellline]["invivo"]+"/chrombpnet_model/chrombpnet_metrics.json"
	tpath = cellline+"/"+data[cellline]["transfer"]+"/chrombpnet_model/chrombpnet_metrics.json"
	idata = json.load(open(ipath))
	tdata = json.load(open(tpath))

	print(cellline)
	invivo.append(idata["counts_metrics"][args.peaks_nonpeaks][args.metric])
	transfer.append(tdata["counts_metrics"][args.peaks_nonpeaks][args.metric])

invivo1 = [np.round(val,2) for val in invivo]
transfer1 = [np.round(val,2) for val in transfer]


print(invivo1)
print(transfer1)

print(stats.ttest_rel(invivo,transfer))


invivo = []
transfer = []

for cellline in data:
	ipath = cellline+"/"+data[cellline]["invivo"]+"/chrombpnet_model/chrombpnet_metrics.json"
	tpath = cellline+"/"+data[cellline]["transfer"]+"/chrombpnet_model/chrombpnet_metrics.json"
	idata = json.load(open(ipath))
	tdata = json.load(open(tpath))

	print(cellline)
	invivo.append(idata["profile_metrics"][args.peaks_nonpeaks]["median_jsd"])
	transfer.append(tdata["profile_metrics"][args.peaks_nonpeaks]["median_jsd"])

invivo1 = [np.round(val,2) for val in invivo]
transfer1 = [np.round(val,2) for val in transfer]


print(invivo1)
print(transfer1)

print(stats.ttest_rel(invivo,transfer))



