#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
import pandas as pd
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_curve
from sklearn.metrics import RocCurveDisplay
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn import metrics
import matplotlib.pyplot as plt
import matplotlib
import random
import os
import argparse
from sklearn.metrics import precision_recall_curve
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve, auc
import json

#data = open("/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/chip-seq/hepg2_chip_match.tsv", "r").readlines()
data = open("/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/chip-seq/gm12878_chip_match.tsv", "r").readlines()
celltype="GM12878"

dictionary = {}
for line in data:

	vals = line.split("\t")

	motif_short=vals[0].strip()
	print(motif_short)

	if len(motif_short) == 0:
		continue
	chip_labels=pd.read_csv("/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/chip-seq/"+motif_short+"/chip_intersect_50bp.bed", sep="\t", header=None)

	print(chip_labels.shape)
	if "FOX" in motif_short:
		motif_name=motif_short+"_"+motif_short+"_1"
	else:
		motif_name=motif_short+"_"+motif_short
	#motif_name="NRF1_HUMAN.H11MO.0.A_NRF1_HUMAN.H11MO.0.A"

	file_id=celltype+"_mean_counts_ATAC_PE"
	moods_dir = "/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/output/"+file_id+"/"+motif_name+"/beds/"
	counts_atac = pd.read_csv(moods_dir+motif_name+"_all.bed", sep="\t", header=None)
	print(counts_atac.shape)
	counts_atac.head()


	file_id=celltype+"_mean_profile_ATAC_PE"
	moods_dir = "/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/output/"+file_id+"/"+motif_name+"/beds/"
	cprofile_atac = pd.read_csv(moods_dir+motif_name+"_all.bed", sep="\t", header=None)


	file_id=celltype+"_mean_counts_DNASE_SE"
	moods_dir = "/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/output/"+file_id+"/"+motif_name+"/beds/"
	counts_dnase = pd.read_csv(moods_dir+motif_name+"_all.bed", sep="\t", header=None)

	file_id=celltype+"_mean_profile_DNASE_SE"
	moods_dir = "/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/output/"+file_id+"/"+motif_name+"/beds/"
	cprofile_dnase = pd.read_csv(moods_dir+motif_name+"_all.bed", sep="\t", header=None)

	file_id=celltype+"_ATAC_TOBIAS"
	moods_dir = "/mnt/lab_data2/anusri/chrombpnet/src/evaluation/figure_notebooks/syntax_stuff/06_25_2022_motif_scanning_tobias_toolkit/"+file_id+"/"+motif_name+"/beds/"
	tobias_atac = pd.read_csv(moods_dir+motif_name+"_all.bed", sep="\t", header=None)

	tobias_atac.head()
	chip_labels.head()

	sum(chip_labels[10]>0)
	chip_labels[11]=chip_labels[10]>0

	counts_dnase.head()

	plt.figure()
	label_col=11
	dictionary[motif_short]={}

	dictionary[motif_short]["TOTAL HITS"] = chip_labels.shape[0]

	dictionary[motif_short]["TOTAL POSITIVE LABELS"] = sum(chip_labels[11])


	fpr_dnase_counts, tpr_dnase_counts, _ = roc_curve(chip_labels[label_col], counts_dnase[9])
	roc_auc = auc(fpr_dnase_counts, tpr_dnase_counts)
	label1="counts dnase AUROC = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["DNASE CHROMBPNET COUNTS"] = np.round(roc_auc,2)
	plt.plot(fpr_dnase_counts, tpr_dnase_counts, label=label1)

	fpr_dnase_profile, tpr_dnase_profile, _ = roc_curve(chip_labels[label_col], cprofile_dnase[9])
	roc_auc = auc(fpr_dnase_profile, tpr_dnase_profile)
	label1="profile dnase AUROC = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["DNASE CHROMBPNET PROFILE"] = np.round(roc_auc,2)

	plt.plot(fpr_dnase_profile, tpr_dnase_profile, label=label1)

	fpr_atac_counts, tpr_atac_counts, _ = roc_curve(chip_labels[label_col], counts_atac[9])
	roc_auc = auc(fpr_atac_counts, tpr_atac_counts)
	label1="counts atac AUROC = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["ATAC CHROMBPNET COUNTS"] = np.round(roc_auc,2)
	plt.plot(fpr_atac_counts, tpr_atac_counts, label=label1)

	fpr_atac_profile, tpr_atac_profile, _ = roc_curve(chip_labels[label_col], cprofile_atac[9])
	roc_auc = auc(fpr_atac_profile, tpr_atac_profile)
	label1="profile atac AUROC = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["ATAC CHROMBPNET PROFILE"] = np.round(roc_auc,2)
	plt.plot(fpr_atac_profile, tpr_atac_profile, label=label1)



	fpr_atac_TOBIAS, tpr_atac_TOBIAS, _ = roc_curve(chip_labels[label_col], tobias_atac[9])
	roc_auc = auc(fpr_atac_TOBIAS, tpr_atac_TOBIAS)
	label1="TOBIAS atac AUROC = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["ATAC TOBIAS"] = np.round(roc_auc,2)

	plt.plot(fpr_atac_TOBIAS, tpr_atac_TOBIAS, label=label1)


	fpr_atac_profile, tpr_atac_profile, _ = roc_curve(chip_labels[label_col], np.max((cprofile_atac[9],counts_atac[9]), axis=0))
	roc_auc = auc(fpr_atac_profile, tpr_atac_profile)
	label1="profile atac AUROC = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["ATAC CHROMBPNET MAX OF PROFILE/COUNTS"] = np.round(roc_auc,2)
	plt.plot(fpr_atac_profile, tpr_atac_profile, label=label1)


	fpr_atac_profile, tpr_atac_profile, _ = roc_curve(chip_labels[label_col], np.max((cprofile_dnase[9],counts_dnase[9]), axis=0))
	roc_auc = auc(fpr_atac_profile, tpr_atac_profile)
	label1="profile atac AUROC = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["DNASE CHROMBPNET MAX OF PROFILE/COUNTS"] = np.round(roc_auc,2)
	plt.plot(fpr_atac_profile, tpr_atac_profile, label=label1)


	fpr_atac_PWM, tpr_atac_PWM, _ = roc_curve(chip_labels[label_col], tobias_atac[4])
	roc_auc = auc(fpr_atac_PWM, tpr_atac_PWM)
	label1="PWM Score = "+str(np.round(roc_auc,2))
	dictionary[motif_short]["PWM"] = np.round(roc_auc,2)
	plt.plot(fpr_atac_PWM, tpr_atac_PWM, label=label1)


	plt.xlabel("False Positive Rate")
	plt.ylabel("True Positive Rate")
	plt.legend()


with open('benchmarking.json', 'w') as fp:
    json.dump(dictionary, fp)
