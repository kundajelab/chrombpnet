import pandas as pd
import os
import numpy as np
from scipy.stats import pearsonr
from scipy.stats import spearmanr
import matplotlib.pyplot as plt

output_path="atac_subsample/"
merged = pd.read_csv(os.path.join(output_path,"variant_scores.tsv"),sep=",", header=0, index_col=False)
merged['logratio'] = merged['logratio'].astype(float)
merged['abs_logratio']=np.abs(merged['logratio'])
merged['abs_Alt_Minus_Ref']=np.abs(merged['Alt_Minus_Ref'])
merged['sig']=merged['pvalue']<5e-5
atac=merged

output_path="dnase_subsample/"
merged = pd.read_csv(os.path.join(output_path,"variant_scores.tsv"),sep=",", header=0, index_col=False)
merged['logratio'] = merged['logratio'].astype(float)
merged['abs_logratio']=np.abs(merged['logratio'])
merged['abs_Alt_Minus_Ref']=np.abs(merged['Alt_Minus_Ref'])
merged['sig']=merged['pvalue']<5e-5
dnase=merged

merged = pd.merge(atac,dnase,on="rsid")
print(merged.head())

#plt.scatter(atac["Alt_Minus_Ref"], dnase["Alt_Minus_Ref"], c="r", alpha=0.5)
plt.scatter(atac["Alt_Minus_Ref"][atac["sig"]==True], dnase["Alt_Minus_Ref"][dnase["sig"]==True], c="r", alpha=0.5)

plt.xlabel("ATAC predictions")
plt.ylabel("DNASE predictions")

plt.savefig("corr_pos.png")
print(pearsonr(merged["Alt_Minus_Ref_x"][merged["sig_x"]==True],merged["Alt_Minus_Ref_y"])[merged["sig_x"]==True])
print(spearmanr(merged["Alt_Minus_Ref_x"][merged["sig_x"]==True],merged["Alt_Minus_Ref_y"][merged["sig_x"]==True]))
