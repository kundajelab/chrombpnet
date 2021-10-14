import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import os
output_path="atac_output"
atac_dsqtl = pd.read_csv(os.path.join(output_path,"variant_scores.tsv"),sep=",", header=0, index_col=False)
output_path="dnase_output"
dnase_dsqtl = pd.read_csv(os.path.join(output_path,"variant_scores.tsv"),sep=",", header=0, index_col=False)


plt.scatter(atac_dsqtl["sum_logratio_pred"][atac_dsqtl["label"]==1], dnase_dsqtl["sum_logratio_pred"][dnase_dsqtl["label"]==1], alpha=0.3, c="r")
#plt.scatter(atac_dsqtl["sum_logratio_pred"], dnase_dsqtl["sum_logratio_pred"], alpha=0.3, c="r")
plt.xlabel("ATAC variant scores")
plt.ylabel("DNASE variant scores")
plt.savefig("correlation_pos.pdf")
#print(pearsonr(atac_dsqtl["sum_logratio_pred"][atac_dsqtl["label"]==1], dnase_dsqtl["sum_logratio_pred"][dnase_dsqtl["label"]==1]))
print(pearsonr(atac_dsqtl["sum_logratio_pred"], dnase_dsqtl["sum_logratio_pred"]))
