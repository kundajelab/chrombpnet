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
data_type="ATAC_PE"
file_id="HEPG2_mean_"+dtype+"_"+data_type
motifs_name="NRF1_HUMAN.H11MO.0.A_NRF1_HUMAN.H11MO.0.A"

moods_dir="/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/output/"+file_id+"/"+motifs_name+"/beds/"
moods_hits = pd.read_csv(moods_dir+motifs_name+"_all.bed", sep="\t", header=None)


bigwig_dir="/mnt/lab_data2/anusri/chrombpnet/results/chrombpnet/"
#bw = pyBigWig.open(bigwig_dir+data_type+"/"+cell_type+"/nautilus_runs_jun16/"+model_name+"/interpret/merged."+cell_type+"."+dtype+".bw")
bw = pyBigWig.open("/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/output/"+file_id.replace("profile","counts")+"/plots/fpd_scores.bw")


scores = []

for i,r in tqdm(moods_hits.iterrows(), total=moods_hits.shape[0]):
    #val = np.nan_to_num(bw.values(r[0],r[1], r[2]))
    #region_score = np.mean(np.abs(val))
    midpoint = (r[1]+r[2])//2
    val = np.nan_to_num(bw.values(r[0],midpoint, midpoint+1))
    region_score = np.mean(val)

    scores.append(region_score)

moods_hits["11"] = scores

moods_hits.to_csv(moods_dir+motifs_name+"_all_with_newvals1.bed", sep="\t", index=False, header=False)





