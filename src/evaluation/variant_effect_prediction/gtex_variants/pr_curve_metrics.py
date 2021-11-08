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

font = {'family' : 'normal',
        'weight' : 'bold',
        'size'   : 11}
matplotlib.rc('font', **font)

parser=argparse.ArgumentParser(description="variant effect scoring on dsQTLS")
parser.add_argument("--output_path")
args = parser.parse_args()

output_path=args.output_path

dsqtl = pd.read_csv(os.path.join(output_path,"variant_scores.tsv"),sep=",", header=0, index_col=False)
dsqtl["pvalue"] = dsqtl["pvalue"].astype(float)
high_sig=dsqtl.loc[dsqtl["pvalue"]<=0.1]
high_sig['label']=True
non_sig=dsqtl.loc[dsqtl["pvalue"]>0.9998]
non_sig['label']=False
dsqtl=pd.concat((high_sig,non_sig),axis=0)


print(dsqtl.head())

fpr_jsd, tpr_jsd, _ = precision_recall_curve(dsqtl["label"], abs(dsqtl["JSD"]))
fpr_jd, tpr_jd, _ = precision_recall_curve(dsqtl["label"], abs(dsqtl["sum_logratio_pred"]))
fpr_gt, tpr_gt, _ = precision_recall_curve(dsqtl["label"], abs(dsqtl["logratio"]))
fpr_counts, tpr_counts, _ = precision_recall_curve(dsqtl["label"], abs(dsqtl["Alt_Minus_Ref"]))

roc_auc = metrics.average_precision_score(fpr_gt, tpr_gt)
plt.plot(fpr_gt,tpr_gt,label="experimental estimate, auc="+str(round(roc_auc,2)))

roc_auc = metrics.average_precision_score(dsqtl["label"], abs(dsqtl["sum_logratio_pred"]))
plt.plot(tpr_jd,fpr_jd,label="sum_logratio difference in profile of ref/alt, auprc="+str(round(roc_auc,2)))

roc_auc = metrics.average_precision_score(dsqtl["label"], abs(dsqtl["JSD"]))
plt.plot(tpr_jsd,fpr_jsd,label="JSD difference in profile of ref/alt, auprc="+str(round(roc_auc,2)))

roc_auc = metrics.average_precision_score(dsqtl["label"], abs(dsqtl["Alt_Minus_Ref"]))
plt.plot(tpr_counts,fpr_counts,label="difference in counts of ref/alt, auprc="+str(round(roc_auc,2)))

#fpr_random, tpr_random, _ = precision_recall_curve(dsqtl["label"], in1)
#roc_auc = metrics.average_precision_score(dsqtl["label"], in1)
#plt.plot(fpr_random,tpr_random,label="random baseline, AP="+str(round(roc_auc,2)))
plt.plot(tpr_gt, [sum(dsqtl["label"]==1)/len(dsqtl["label"])]*len(tpr_gt), linewidth=2, label="random baseline, AP="+str(np.round(sum(dsqtl["label"]==1)/len(dsqtl["label"]),2)))

plt.legend(loc='upper right')
plt.xlabel("Recall")
plt.ylabel("Precision")
plt.savefig(os.path.join(output_path,"pr_curve_metrics.pdf"))

