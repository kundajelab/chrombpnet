rm(list=ls())
library(ggplot2)
library(ggridges)
library(coin)

#tf='spi1'
tf='junD'
data=read.table(paste("5e-5.gm12878.",tf,".tsv",sep=''),header=TRUE,sep=',')
data$sig=data$pvalue<5e-5

data$sig=factor(data$sig)
sig_logratio=data$sum_logratio_pred[data$sig==TRUE]
num_sig=length(sig_logratio)
non_sig_logratio=data$sum_logratio_pred[data$sig==FALSE]
num_non_sig=length(non_sig_logratio)
#ks.test(sig_delta,non_sig_delta)
wilcox_test(sum_logratio_pred~sig,data=data,conf.level=0.99,conf.int=TRUE)
ggplot(data, aes(sum_logratio_pred, colour =sig)) +
  stat_ecdf(size=1)+
  scale_color_manual(values=c("#e41a1c","#377eb8a0"))+
  xlab("Abs(Log(Prob(Alt))-Log(Prob(Ref)))")+
  ylab("CDF")+
  ggtitle(paste(tf, ' p-val<5e-5 (n=',num_sig,' \n vs p-val>5e-5 (n=',num_non_sig,sep=''))+
  theme_bw(15)+
  theme(legend.position="top")
