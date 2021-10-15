import os
import argparse
from plotnine import * 
import scipy
import pandas as pd
import numpy as np

parser=argparse.ArgumentParser(description="variant effect scoring on dsQTLS")
parser.add_argument("--output_path")
args = parser.parse_args()

output_path=args.output_path

merged = pd.read_csv(os.path.join(output_path,"variant_scores.tsv"),sep=",", header=0, index_col=False)
merged['logratio'] = merged['logratio'].astype(float)
merged['abs_logratio']=np.abs(merged['logratio'])
merged['abs_Alt_Minus_Ref']=np.abs(merged['Alt_Minus_Ref'])
merged['sig']=merged['pvalue']<5e-5

print(scipy.stats.pearsonr(merged['abs_logratio'],merged['JSD']))
print(scipy.stats.pearsonr(merged['abs_logratio'],merged['sum_logratio_pred']))
print(scipy.stats.pearsonr(merged['logratio'],merged['Alt_Minus_Ref']))

print(scipy.stats.spearmanr(merged['abs_logratio'],merged['JSD']))
print(scipy.stats.spearmanr(merged['abs_logratio'],merged['sum_logratio_pred']))
print(scipy.stats.spearmanr(merged['logratio'],merged['Alt_Minus_Ref']))

p = (ggplot(merged, aes('logratio', 'Alt_Minus_Ref', color='factor(sig)'))
 + geom_point(alpha=1.0,size=1.0)
 + xlab("logratio bQTL")
 + ylab("Log(Counts(Alt)) -\nLog(Counts(Ref)))")
 + geom_smooth(aes(group=1),colour="black",method="lm")
 + scale_color_manual(values=['#e41a1c','#377eb8'],name="Significant bQTL?")
 + theme_bw(20))

p.save(os.path.join(output_path,"scatter_top.png"))
