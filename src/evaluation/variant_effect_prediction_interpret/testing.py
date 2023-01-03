import pandas as pd
from plotnine import *
import scipy

snp_regions=pd.read_csv("variant_scores.tsv",header=0,sep='\t')
snp_regions[["logratio", "sig"]] = snp_regions["META_DATA"].str.split(',', expand=True).astype(float)

print(scipy.stats.pearsonr(snp_regions['logratio'],snp_regions['log_probs_diff_abs_sum']))

p = (ggplot(snp_regions, aes('logratio', 'log_probs_diff_abs_sum', color='factor(sig)'))
 + geom_point(alpha=1.0,size=1.0)
 + xlab("logratio bQTL")
 + ylab("Log(Counts(Alt)) -\nLog(Counts(Ref)))")
 + geom_smooth(aes(group=1),colour="black",method="lm")
 + scale_color_manual(values=['#e41a1c','#377eb8'],name="Significant bQTL?")
 + theme_bw(20))

print()
p.save("scatter_top.png")

