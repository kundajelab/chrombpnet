import deepdish
from statsmodels.distributions.empirical_distribution import ECDF
import matplotlib.pyplot as plt
import pandas as pd
import pyBigWig
import numpy as np
from tqdm import tqdm


cell_type="HEPG2"
model_name="HEPG2_05.09.2022_bias_128_4_1234_0.8_fold_0"
dtype="profile"
#data_type="DNASE_PE"
data_type="ATAC_PE"
dir="/mnt/lab_data3/anusri/chrombpnet/results/chrombpnet/"+data_type+"/"+cell_type+"/"+model_name+"/06_22_2022_motif_scanning/"
bigwig_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"
bw = pyBigWig.open(bigwig_dir+data_type+"/"+cell_type+"/nautilus_runs_jun16/"+model_name+"/interpret/merged."+cell_type+"."+dtype+".bw")

hdf5 = dir+"/interpret/5K.negatives."+dtype+"_scores.h5"
moods_hits = pd.read_csv(dir+"/moods/"+dtype+"/moods_scan_on_modisco_full.bed", sep="\t", header=None)


print(moods_hits.head())

d = deepdish.io.load(hdf5, '/projected_shap/seq')
deepshaps_in_neg = d.flatten()

print("quantile 0.95: ",np.quantile(deepshaps_in_neg, 0.95))
print("quantile 0.05: ", np.quantile(deepshaps_in_neg, 0.05))

ecdf = ECDF(deepshaps_in_neg)

scores = []
pvalues = []

for i,r in tqdm(moods_hits.iterrows(), total=moods_hits.shape[0]):
    val = np.nan_to_num(bw.values(r[0],r[1], r[2]))
    region_score = np.mean(val)
    scores.append(region_score)
    pvalues.append(ecdf(region_score))

moods_hits["11"] = scores
moods_hits["12"] = pvalues

moods_hits.to_csv(dir+"/moods/"+dtype+"/moods_modisco_hits_with_pvalues.bed", sep="\t", index=False, header=False)




