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
print(dsqtl.head())

fpr_chrombpnet, tpr_chrombpnet, _ = precision_recall_curve(dsqtl["label"], abs(dsqtl["sum_logratio_pred"]))
roc_auc = metrics.average_precision_score(dsqtl["label"], abs(dsqtl["sum_logratio_pred"]))
plt.plot(tpr_chrombpnet,fpr_chrombpnet,  linewidth=2, label="Chrom-BPnet, AP="+str(round(roc_auc,2)))

#fpr_svm, tpr_svm, _ = precision_recall_curve(dsqtl["label"], abs(dsqtl["svm_pred"]))
#roc_auc = metrics.average_precision_score(dsqtl["label"], abs(dsqtl["svm_pred"]))
#plt.plot(tpr_svm, fpr_svm, label="center-weighted gkm SVM, AP="+str(round(roc_auc,2)))

fpr_gkm, tpr_gkm, _ = precision_recall_curve(dsqtl["label"], dsqtl["abs_gkm_SVM"])
roc_auc = metrics.average_precision_score(dsqtl["label"], dsqtl["abs_gkm_SVM"])
#plt.plot(tpr_gkm, fpr_gkm, linewidth=2, label="gkm SVM (from deltaSVM paper), AP="+str(round(roc_auc,2)))
plt.plot(tpr_gkm, fpr_gkm, linewidth=2, label="gkm SVM (from deltaSVM paper), AP="+str(round(roc_auc,2)))

in1 = list(dsqtl["label"])
random.shuffle(in1)
in1=np.array(in1)

fpr_random, tpr_random, _ = precision_recall_curve(dsqtl["label"], in1)
roc_auc = metrics.average_precision_score(dsqtl["label"], in1)
#plt.plot(fpr_random,tpr_random,label="random baseline, AP="+str(round(roc_auc,2)))
plt.plot(tpr_random, [sum(dsqtl["label"]==1)/len(dsqtl["label"])]*len(fpr_random), linewidth=2, label="random baseline, AP="+str(round(roc_auc,2)))

print(len(dsqtl["label"]), sum(dsqtl["label"]==1))
#print(dsqtl["label"])
plt.legend(loc='upper right')
plt.xlabel("Recall")
plt.ylabel("Precision")
#pplt.show()
plt.savefig(os.path.join(output_path,"pr_curve.pdf"))

